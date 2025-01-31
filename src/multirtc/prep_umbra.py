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
    def calculate_instantaneous_prf(time: float, ipp):
        if not ipp.TStart <= time <= ipp.TEnd:
            raise ValueError('Time must be within the IPP')
        check_poly_order(ipp.IPPPoly)
        prf_coeff = np.polyder(ipp.IPPPoly.Coefs[::-1])
        prf = np.polyval(prf_coeff, time)
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
        scp_tcoa = sicd.Grid.TimeCOAPoly.Coefs[0, 0]
        ipp = list(sicd.Timeline.IPP)[0]
        prf = cls.calculate_instantaneous_prf(scp_tcoa, ipp)
        orbit = cls.calculate_orbit(
            sensing_period_start,
            sicd.ImageFormation.TStartProc,
            scp_tcoa,
            sicd.ImageFormation.TEndProc,
            sicd.Position.ARPPoly,
        )
        umbra_sicd = cls(
            id=sicd.CollectionInfo.CoreName,
            file_path=file_path,
            wavelength=wavelength,
            lookside=lookside,
            starting_range=sicd.SCPCOA.SlantRange,  # Range at scene center, not starting range
            prf=prf,
            # prf=486.0,
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
        )
        return umbra_sicd

    def as_isce3_radargrid(self):
        radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start_sec,
            wavelength=self.wavelength,
            lookside=self.lookside,
            length=self.shape[0],
            width=self.shape[1],
            ref_epoch=isce3.core.DateTime(self.reference_epoch),
            range_pixel_spacing=self.range_step,
            starting_range=self.starting_range,
            prf=self.prf,
            # range_pixel_spacing=2.3,
            # prf=486,
            # starting_range=974673,
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
