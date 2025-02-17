import argparse
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

import isce3
import numpy as np

# from numpy.polynomial.polynomial import polyval
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Polygon

from multirtc import dem


def check_poly_order(poly):
    assert len(poly.Coefs) == poly.order1 + 1, 'Polynomial order does not match number of coefficients'


@dataclass
class Point:
    x: float
    y: float


@dataclass
class UmbraSICD:
    id: str
    file_path: Path
    sicd: str
    wavelength: float
    lookside: str  # 'right' or 'left'
    starting_range: float
    end_range: float
    prf: float
    range_step: float  # meters
    azimuth_step: float  # meters
    beta0_coeff: np.ndarray
    orbit: isce3.core.Orbit
    row_uvect: np.ndarray
    col_uvect: np.ndarray
    arp_pos: np.ndarray
    arp_vel: np.ndarray
    reference_epoch: datetime
    sensing_start: datetime
    sensing_start_sec: float
    sensing_end_sec: float
    shape: tuple[int, int]  # (rows, cols)
    scp_index: tuple[int, int]  # (rows, cols)
    scp_pos: np.ndarray
    footprint: Polygon
    center: Point

    def load_data(self):
        """Load data from the UMBRA SICD file."""
        reader = SICDReader(str(self.file_path))
        data = reader[:, :]
        return data

    def get_xrow_ycol(self) -> tuple[np.ndarray]:
        """Calculate xrow and ycol for the umbra_sicd."""
        shadows_down_n_rows = self.shape[0]
        shadows_down_n_cols = self.shape[1]

        irow = np.tile(np.arange(shadows_down_n_cols), (shadows_down_n_rows, 1)).T
        irow -= self.scp_index[1]
        xrow = irow * self.range_step

        icol = np.tile(np.arange(shadows_down_n_rows), (shadows_down_n_cols, 1))
        icol -= self.scp_index[0]
        ycol = icol * self.azimuth_step
        return xrow, ycol

    def get_doppler_centroid_grid(self, n_samples=50, pixel_buffer=1_000):
        half_span = (self.shape[0] // 2) + pixel_buffer
        rows = np.linspace(-half_span, half_span, n_samples, dtype=int)
        row_offset = ((rows * self.sicd.Grid.Row.SS)[:, np.newaxis] * self.row_uvect).T
        row_ecef = row_offset + self.scp_pos[:, np.newaxis]
        row_vec = row_ecef - self.arp_pos[:, np.newaxis]
        # row_vec[1, :] -= 34_000
        row_mag = np.linalg.norm(row_vec, axis=0)
        row_look = row_vec / np.linalg.norm(row_vec, axis=0)

        v_mag = np.linalg.norm(self.arp_vel)
        v_hat = self.arp_vel / v_mag
        row_squint = np.arcsin(np.dot(row_look.T, v_hat))
        row_dc = 2.0 / self.wavelength * v_mag * np.sin(row_squint)

        doppler = np.tile(row_dc, (n_samples, 1))
        avg_diff = np.mean(np.diff(row_mag))
        ranges = np.arange(row_mag[0], row_mag[-1] + 0.01, avg_diff)
        azimuth_times = np.linspace(-0.75, 2.75, n_samples) * self.sicd.SCPCOA.SCPTime

        doppler_lut = isce3.core.LUT2d(xcoord=ranges, ycoord=azimuth_times, data=doppler)
        return doppler_lut

    @staticmethod
    def calculate_orbit(sensing_period_start, position):
        pos_poly = [position.ARPPoly.X, position.ARPPoly.Y, position.ARPPoly.Z]
        [check_poly_order(poly) for poly in pos_poly]
        pos_poly = [poly[::-1] for poly in pos_poly]
        vel_poly = [np.polyder(poly) for poly in pos_poly]
        times = np.arange(-30, 31, 1)
        svs = []
        for time in times:
            pos = [np.polyval(poly, time) for poly in pos_poly]
            vel = [np.polyval(poly, time) for poly in vel_poly]
            dt_time = isce3.core.DateTime(sensing_period_start + timedelta(seconds=float(time)))
            sv = isce3.core.StateVector(dt_time, pos, vel)
            svs.append(sv)
        orbit = isce3.core.Orbit(svs, isce3.core.DateTime(sensing_period_start))
        orbit.set_interp_method('Legendre')
        return orbit

    @staticmethod
    def calculate_instantaneous_prf(time: float, ipp, sicd):
        if not ipp.TStart <= time <= ipp.TEnd:
            raise ValueError('Time must be within the IPP')
        check_poly_order(ipp.IPPPoly)
        prf_coeff = np.polyder(ipp.IPPPoly.Coefs[::-1])
        prf = np.polyval(prf_coeff, time)
        return prf

    @staticmethod
    def calculate_start_end_range(sicd):
        sicd_row_uvect = sicd.Grid.Row.UVectECF
        row_uvect = np.array([sicd_row_uvect.X, sicd_row_uvect.Y, sicd_row_uvect.Z])
        row_offset = sicd.ImageData.SCPPixel.Row * sicd.Grid.Row.SS * row_uvect

        scp_ecf = np.array([sicd.GeoData.SCP.ECF.X, sicd.GeoData.SCP.ECF.Y, sicd.GeoData.SCP.ECF.Z])
        arp_ecf = np.array([sicd.SCPCOA.ARPPos.X, sicd.SCPCOA.ARPPos.Y, sicd.SCPCOA.ARPPos.Z])

        start_range_vec = scp_ecf - row_offset - arp_ecf
        start_range_mag = np.linalg.norm(start_range_vec, axis=0)

        end_range_vec = scp_ecf - row_offset - arp_ecf
        end_range_mag = np.linalg.norm(end_range_vec, axis=0)
        return start_range_mag, end_range_mag

    def get_stripmap_prf(self):
        vel_along_track = np.dot(self.arp_vel, self.col_uvect)
        vel_mag = np.linalg.norm(vel_along_track)
        prf = vel_mag / self.sicd.Grid.Col.SS
        return prf

    @classmethod
    def from_sarpy_sicd(cls, sicd, file_path):
        center_frequency = sicd.RadarCollection.TxFrequency.Min + sicd.RadarCollection.TxFrequency.Max / 2
        wavelength = isce3.core.speed_of_light / center_frequency
        lookside = isce3.core.LookSide.Right if sicd.SCPCOA.SideOfTrack == 'R' else isce3.core.LookSide.Left
        sensing_period_start = sicd.Timeline.CollectStart.astype('M8[ms]').astype('O')
        range_step = sicd.Grid.Row.SS
        azimuth_step = sicd.Grid.Col.SS
        footprint = Polygon([(ic.Lon, ic.Lat) for ic in sicd.GeoData.ImageCorners])
        scp_coa = sicd.SCPCOA
        ipp = list(sicd.Timeline.IPP)[0]
        prf = cls.calculate_instantaneous_prf(scp_coa.SCPTime, ipp, sicd)
        orbit = cls.calculate_orbit(sensing_period_start, sicd.Position)
        starting_range, end_range = cls.calculate_start_end_range(sicd)
        umbra_sicd = cls(
            id=sicd.CollectionInfo.CoreName,
            file_path=file_path,
            sicd=sicd,
            wavelength=wavelength,
            lookside=lookside,
            starting_range=starting_range,
            prf=prf,
            end_range=end_range,
            range_step=range_step,
            azimuth_step=azimuth_step,
            beta0_coeff=sicd.Radiometric.BetaZeroSFPoly.Coefs,
            orbit=orbit,
            reference_epoch=sensing_period_start,
            sensing_start=sensing_period_start + timedelta(seconds=sicd.ImageFormation.TStartProc),
            sensing_start_sec=sicd.ImageFormation.TStartProc,
            sensing_end_sec=sicd.ImageFormation.TEndProc,
            row_uvect=np.array([sicd.Grid.Row.UVectECF.X, sicd.Grid.Row.UVectECF.Y, sicd.Grid.Row.UVectECF.Z]),
            col_uvect=np.array([sicd.Grid.Col.UVectECF.X, sicd.Grid.Col.UVectECF.Y, sicd.Grid.Col.UVectECF.Z]),
            arp_vel=np.array([sicd.SCPCOA.ARPVel.X, sicd.SCPCOA.ARPVel.Y, sicd.SCPCOA.ARPVel.Z]),
            arp_pos=np.array([sicd.SCPCOA.ARPPos.X, sicd.SCPCOA.ARPPos.Y, sicd.SCPCOA.ARPPos.Z]),
            shape=(sicd.ImageData.NumRows, sicd.ImageData.NumCols),  # backwards for shadows-down
            scp_index=(sicd.ImageData.SCPPixel.Row, sicd.ImageData.SCPPixel.Col),  # backwards for shadows-down
            scp_pos=np.array([sicd.GeoData.SCP.ECF.X, sicd.GeoData.SCP.ECF.Y, sicd.GeoData.SCP.ECF.Z]),
            footprint=footprint,
            center=Point(sicd.GeoData.SCP.LLH.Lat, sicd.GeoData.SCP.LLH.Lon),
        )
        return umbra_sicd

    def as_isce3_radargrid(self):
        radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start_sec,
            wavelength=self.wavelength,
            lookside=self.lookside,
            # length/width are backwards for shadows-down
            length=self.shape[1],
            width=self.shape[0],
            ref_epoch=isce3.core.DateTime(self.reference_epoch),
            range_pixel_spacing=self.range_step,
            starting_range=self.starting_range,
            prf=self.get_stripmap_prf(),
        )
        assert radar_grid.ref_epoch == self.orbit.reference_epoch
        return radar_grid


def prep_umbra(granule_path: Path, work_dir: Optional[Path] = None) -> Path:
    """Prepare data for burst-based processing.

    Args:
        granule_path: Path to the UMBRA SICD file
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    reader = SICDReader(str(granule_path))
    sicd = reader.get_sicds_as_tuple()[0]
    umbra_sicd = UmbraSICD.from_sarpy_sicd(sicd, granule_path)

    dem_path = work_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, umbra_sicd.footprint)
    return umbra_sicd, dem_path


def main():
    """Prep SLC entrypoint.

    Example command:
    prep_burst CR-28_2024-12-03-18-21-21_UMBRA-10_SICD_MM.nitf
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('granule', help='Umbra SICD to load data for.')
    parser.add_argument('--work-dir', default=None, help='Working directory for processing')

    args = parser.parse_args()

    prep_umbra(**args.__dict__)


if __name__ == '__main__':
    main()
