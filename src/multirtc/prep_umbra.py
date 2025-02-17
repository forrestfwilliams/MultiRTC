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
    wavelength: float
    lookside: str  # 'right' or 'left'
    starting_range: float
    prf: float
    range_step: float  # meters
    azimuth_step: float  # meters
    beta0_coeff: np.ndarray
    orbit: isce3.core.Orbit
    reference_epoch: datetime
    sensing_start: datetime
    sensing_start_sec: float
    sensing_end_sec: float
    shape: tuple[int, int]  # (rows, cols)
    scp_index: tuple[int, int]  # (rows, cols)
    footprint: Polygon
    center: Point
    arp_pos: np.ndarray
    arp_vel: np.ndarray
    row_uvect: np.ndarray
    row_ss: float
    col_uvect: np.ndarray
    col_ss: float

    def load_data(self):
        """Load data from the UMBRA SICD file."""
        reader = SICDReader(str(self.file_path))
        data = reader[:, :]
        return data

    @staticmethod
    def calculate_orbit(
        sensing_period_start: datetime, sensing_start: float, scp_tcoa: float, sensing_end: float, state_poly
    ):
        """Calculate the orbit for a sicd.
        isce3.core.Orbit takes two arguments, a set of orbit state vectors and a reference epoch.

        The reference epoch is starting time - 2 days.
        State vectors is a list of isce3.core.StateVector objects where each object takes the form:
        isce3.core.StateVector(isce3_datetime, [x_pos, y_pos, z_pos], [x_vel, y_vel, z_vel])

        Spotlight images have a constant azimuth time within a range line, so only one state vector is needed???
        """
        # FIXME: Use the sicd.SCPCOA instead
        pos_scp, vel_scp = [], []
        for poly in [state_poly.X, state_poly.Y, state_poly.Z]:
            check_poly_order(poly)
            pos_scp.append(np.polyval(poly.Coefs[::-1], scp_tcoa))
            vel_coeff = np.polyder(poly.Coefs[::-1])
            vel_scp.append(np.polyval(vel_coeff, scp_tcoa))

        # Umbra only gives us one state vector, so assume a constant velocity to get the rest of the local orbit
        total_time = sensing_end - sensing_start
        half_time = total_time // 2
        times_relative = np.arange(-half_time, half_time, 1)
        times = times_relative + scp_tcoa
        svs = []
        for time in times:
            pos = list(np.array(pos_scp) + (np.array(vel_scp) * time))
            sv_time = isce3.core.DateTime(sensing_period_start + timedelta(seconds=time))
            svs.append(isce3.core.StateVector(sv_time, pos, vel_scp))

        orbit = isce3.core.Orbit(svs, isce3.core.DateTime(sensing_period_start))
        return orbit

    @staticmethod
    def calculate_orbit_v2(scp_coa, sensing_period_start):
        time = np.arange(-10, 11, 1)

        x_vel = scp_coa.ARPVel.X + scp_coa.ARPAcc.X * time
        x_pos = scp_coa.ARPPos.X + scp_coa.ARPVel.X * time + 0.5 * scp_coa.ARPAcc.X * time**2

        y_vel = scp_coa.ARPVel.Y + scp_coa.ARPAcc.Y * time
        y_pos = scp_coa.ARPPos.Y + scp_coa.ARPVel.Y * time + 0.5 * scp_coa.ARPAcc.Y * time**2

        z_vel = scp_coa.ARPVel.Z + scp_coa.ARPAcc.Z * time
        z_pos = scp_coa.ARPPos.Z + scp_coa.ARPVel.Z * time + 0.5 * scp_coa.ARPAcc.Z * time**2

        time_sensing_period = time + scp_coa.SCPTime

        svs = []
        for i in range(len(time)):
            sv_time = isce3.core.DateTime(sensing_period_start + timedelta(seconds=time_sensing_period[i]))
            sv = isce3.core.StateVector(sv_time, [x_pos[i], y_pos[i], z_pos[i]], [x_vel[i], y_vel[i], z_vel[i]])
            svs.append(sv)

        orbit = isce3.core.Orbit(svs, isce3.core.DateTime(sensing_period_start))
        orbit.set_interp_method('Legendre')
        return orbit

    @staticmethod
    def calculate_instantaneous_prf(time: float, ipp):
        if not ipp.TStart <= time <= ipp.TEnd:
            raise ValueError('Time must be within the IPP')
        check_poly_order(ipp.IPPPoly)
        prf_coeff = np.polyder(ipp.IPPPoly.Coefs[::-1])
        prf = np.polyval(prf_coeff, time)
        return prf

    @staticmethod
    def calculate_starting_range(sicd):
        scp_range = sicd.SCPCOA.SlantRange
        grazing_angle = np.deg2rad(sicd.SCPCOA.GrazeAng)  # TODO: is this for slant range?
        starting_range_to_scene_center = sicd.Grid.Col.SS * sicd.ImageData.SCPPixel.Col
        starting_range = np.sqrt(
            scp_range**2
            + starting_range_to_scene_center**2
            - (2 * scp_range * starting_range_to_scene_center * np.cos(grazing_angle))
        )
        return starting_range

    def get_stripmap_prf(self):
        vel_along_track = np.dot(self.arp_vel, self.col_uvect)
        vel_mag = np.linalg.norm(vel_along_track)
        prf = vel_mag / self.col_ss
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
        prf = cls.calculate_instantaneous_prf(scp_coa.SCPTime, ipp)
        orbit = cls.calculate_orbit_v2(scp_coa, sensing_period_start)
        starting_range = cls.calculate_starting_range(sicd)
        umbra_sicd = cls(
            id=sicd.CollectionInfo.CoreName,
            file_path=file_path,
            wavelength=wavelength,
            lookside=lookside,
            starting_range=starting_range,
            prf=prf,
            range_step=range_step,
            azimuth_step=azimuth_step,
            beta0_coeff=sicd.Radiometric.BetaZeroSFPoly.Coefs,
            orbit=orbit,
            reference_epoch=sensing_period_start,
            sensing_start=sensing_period_start + timedelta(seconds=sicd.ImageFormation.TStartProc),
            sensing_start_sec=sicd.ImageFormation.TStartProc,
            sensing_end_sec=sicd.ImageFormation.TEndProc,
            shape=(sicd.ImageData.NumRows, sicd.ImageData.NumCols),
            scp_index=(sicd.ImageData.SCPPixel.Row, sicd.ImageData.SCPPixel.Col),
            footprint=footprint,
            center=Point(sicd.GeoData.SCP.LLH.Lon, sicd.GeoData.SCP.LLH.Lat),
            arp_pos=np.array([sicd.SCPCOA.ARPPos.X, sicd.SCPCOA.ARPPos.Y, sicd.SCPCOA.ARPPos.Z]),
            arp_vel=np.array([sicd.SCPCOA.ARPVel.X, sicd.SCPCOA.ARPVel.Y, sicd.SCPCOA.ARPVel.Z]),
            row_uvect=np.array([sicd.Grid.Row.UVectECF.X, sicd.Grid.Row.UVectECF.Y, sicd.Grid.Row.UVectECF.Z]),
            row_ss=sicd.Grid.Row.SS,
            col_uvect=np.array([sicd.Grid.Col.UVectECF.X, sicd.Grid.Col.UVectECF.Y, sicd.Grid.Col.UVectECF.Z]),
            col_ss=sicd.Grid.Col.SS,
        )
        return umbra_sicd

    def as_isce3_radargrid(self):
        radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start_sec,
            wavelength=self.wavelength,
            lookside=self.lookside,
            length=self.shape[1],
            width=self.shape[0],
            ref_epoch=isce3.core.DateTime(self.reference_epoch),
            range_pixel_spacing=self.range_step,
            starting_range=self.starting_range,
            prf=self.get_stripmap_prf()
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
