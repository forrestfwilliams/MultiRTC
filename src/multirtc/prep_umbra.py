import argparse
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

import isce3
import numpy as np
from numpy.polynomial.polynomial import polyval
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Polygon

from multirtc import dem


@dataclass
class UmbraSICD:
    id: str
    wavelength: float
    lookside: str  # 'right' or 'left'
    starting_range: float
    prf: float
    range_step: float
    beta0_coeff: np.ndarray
    orbit: isce3.core.Orbit
    sensing_start: datetime
    sensing_end: datetime
    shape: tuple[int, int]  # (rows, cols)
    footprint: Polygon

    @staticmethod
    def calculate_orbit(sensing_start: datetime, sensing_end: datetime, scp_tcoa: int, state_poly):
        """Calculate the orbit for a sicd.
        isce3.core.Orbit takes two arguments, a set of orbit state vectors and a reference epoch.

        The reference epoch is starting time - 2 days.
        State vectors is a list of isce3.core.StateVector objects where each object takes the form:
        isce3.core.StateVector(isce3_datetime, [x_pos, y_pos, z_pos], [x_vel, y_vel, z_vel])

        Spotlight images have a constant azimuth time within a range line, so only one state vector is needed???
        """
        pos, vel = [], []
        for poly in [state_poly.X, state_poly.Y, state_poly.Z]:
            assert len(poly.Coefs) == poly.order1 + 1, 'Polynomial order does not match number of coefficients'
            pos.append(polyval(scp_tcoa, poly.Coefs))
            vel_coeff = np.polyder(poly.Coefs[::-1])[::-1]
            vel.append(polyval(scp_tcoa, vel_coeff))

        sv_start = isce3.core.StateVector(isce3.core.DateTime(sensing_start - timedelta(minutes=1)), pos, vel)
        sv_end = isce3.core.StateVector(isce3.core.DateTime(sensing_end + timedelta(minutes=1)), pos, vel)
        return isce3.core.Orbit([sv_start, sv_end], isce3.core.DateTime(sensing_start - timedelta(minutes=5)))

    @classmethod
    def from_sarpy_sicd(cls, sicd):
        center_frequency = sicd.RadarCollection.TxFrequency.Min + sicd.RadarCollection.TxFrequency.Max / 2
        wavelength = isce3.core.speed_of_light / center_frequency
        lookside = isce3.core.LookSide.Right if sicd.SCPCOA.SideOfTrack == 'R' else isce3.core.LookSide.Left
        sensing_period_start = sicd.Timeline.CollectStart.astype('M8[ms]').astype('O')
        sensing_start = sensing_period_start + timedelta(sicd.ImageFormation.TStartProc)
        sensing_end = sensing_period_start + timedelta(sicd.ImageFormation.TEndProc)
        ipp = list(sicd.Timeline.IPP)[0]
        prf = (ipp.IPPEnd - ipp.IPPStart) / (ipp.TEnd - ipp.TStart)  # not sure if this is correct
        range_step = sicd.Grid.Row.SS
        footprint = Polygon([(ic.Lon, ic.Lat) for ic in sicd.GeoData.ImageCorners])
        scp_tcoa = sicd.Grid.TimeCOAPoly.Coefs[0, 0]
        orbit = cls.calculate_orbit(sensing_start, sensing_end, scp_tcoa, sicd.Position.ARPPoly)
        umbra_sicd = cls(
            id=sicd.CollectionInfo.CoreName,
            wavelength=wavelength,
            lookside=lookside,
            starting_range=sicd.SCPCOA.SlantRange,  # Range at scene center, not starting range
            prf=prf,
            range_step=range_step,
            beta0_coeff=sicd.Radiometric.BetaZeroSFPoly.Coefs,
            orbit=orbit,
            sensing_start=isce3.core.DateTime(sensing_start),
            sensing_end=isce3.core.DateTime(sensing_end),
            shape=(sicd.ImageData.NumRows, sicd.ImageData.NumCols),
            footprint=footprint,
        )
        return umbra_sicd


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
    umbra_sicd = UmbraSICD.from_sarpy_sicd(sicd)

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
