from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

import isce3
import numpy as np
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Polygon

from multirtc import dem


def check_poly_order(poly):
    assert len(poly.Coefs) == poly.order1 + 1, 'Polynomial order does not match number of coefficients'


def to_isce_datetime(dt):
    if isinstance(dt, datetime):
        return isce3.core.DateTime(dt)
    elif isinstance(dt, np.datetime64):
        return isce3.core.DateTime(dt.item())
    else:
        raise ValueError(f'Unsupported datetime type: {type(dt)}. Expected datetime or np.datetime64.')


@dataclass
class Point:
    x: float
    y: float


@dataclass
class CapellaSICD:
    id: str
    file_path: Path
    footprint: Polygon
    wavelength: float
    polarization: str
    lookside: str  # 'right' or 'left'
    hae: float
    shape: tuple[int, int]
    scp_index: tuple[int, int]
    range_pixel_spacing: float
    collect_start: datetime
    sensing_start: float
    sensing_end: float
    prf: float
    starting_range: float
    orbit: isce3.core.Orbit
    beta0_coeff: np.ndarray

    @staticmethod
    def calculate_orbit(epoch, sensing_start, sensing_end, arp_pos_poly):
        svs = []
        orbit_start = np.floor(sensing_start) - 5
        orbit_end = np.ceil(sensing_end) + 5
        for offset_sec in np.arange(orbit_start, orbit_end + 1, 1):
            t = sensing_start + offset_sec
            pos = arp_pos_poly(t)
            vel = arp_pos_poly.derivative_eval(t)
            t_isce = to_isce_datetime(epoch + timedelta(seconds=t))
            svs.append(isce3.core.StateVector(t_isce, pos, vel))
        return isce3.core.Orbit(svs, to_isce_datetime(epoch))

    def as_isce3_radar_grid(self):
        radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start,
            wavelength=self.wavelength,
            prf=self.prf,
            starting_range=self.starting_range,
            range_pixel_spacing=self.range_pixel_spacing,
            lookside=isce3.core.LookSide.Right if self.lookside == 'right' else isce3.core.LookSide.Left,
            length=self.shape[1],  # flipped for "shadows down" convention
            width=self.shape[0],  # flipped for "shadows down" convention
            ref_epoch=to_isce_datetime(self.collect_start),
        )
        return radar_grid

    def get_doppler_centroid_grid(self):
        return isce3.core.LUT2d()

    @classmethod
    def from_sarpy_sicd(cls, sicd, file_path):
        center_frequency = sicd.RadarCollection.TxFrequency.Min + sicd.RadarCollection.TxFrequency.Max / 2
        wavelength = isce3.core.speed_of_light / center_frequency
        polarization = sicd.RadarCollection.RcvChannels[0].TxRcvPolarization.replace(':', '')
        lookside = 'right' if sicd.SCPCOA.SideOfTrack == 'R' else 'left'
        footprint = Polygon([(ic.Lon, ic.Lat) for ic in sicd.GeoData.ImageCorners])
        shape = np.array([sicd.ImageData.NumRows, sicd.ImageData.NumCols])
        scp_index = np.array([sicd.ImageData.SCPPixel.Row, sicd.ImageData.SCPPixel.Col])
        row_shift = sicd.ImageData.SCPPixel.Row - sicd.ImageData.FirstRow
        col_shift = sicd.ImageData.SCPPixel.Col - sicd.ImageData.FirstCol
        row_mult = sicd.Grid.Row.SS
        col_mult = sicd.Grid.Col.SS
        range_pixel_spacing = row_mult
        assert sicd.Grid.Type == 'RGZERO', 'Only range zero doppler grids supported for Capella data'
        collect_start = sicd.Timeline.CollectStart
        # seconds after collect start
        first_col_time = sicd.Grid.TimeCOAPoly((0 - row_shift) * row_mult, (0 - col_shift) * col_mult)
        last_col_time = sicd.Grid.TimeCOAPoly((0 - row_shift) * row_mult, (shape[1] - col_shift) * col_mult)
        sensing_start = min(first_col_time, last_col_time)  # + 2  # fudged
        sensing_end = max(first_col_time, last_col_time)
        prf = np.mean([ipp.IPPPoly.derivative_eval((ipp.TStart + ipp.TEnd) / 2) for ipp in sicd.Timeline.IPP])
        starting_row_pos = (
            sicd.GeoData.SCP.ECF.get_array() + sicd.Grid.Row.UVectECF.get_array() * (0 - row_shift) * row_mult
        )
        starting_range = np.linalg.norm(sicd.SCPCOA.ARPPos.get_array() - starting_row_pos)
        orbit = CapellaSICD.calculate_orbit(collect_start, sensing_start, sensing_end, sicd.Position.ARPPoly)
        beta0_coeff = sicd.Radiometric.BetaZeroSFPoly.Coefs
        capella_sicd = cls(
            id=Path(file_path).with_suffix('').name,
            file_path=file_path,
            footprint=footprint,
            shape=shape,
            wavelength=wavelength,
            polarization=polarization,
            lookside=lookside,
            hae=float(sicd.GeoData.SCP.LLH.HAE),
            range_pixel_spacing=range_pixel_spacing,
            scp_index=scp_index,
            collect_start=collect_start,
            sensing_start=sensing_start,
            sensing_end=sensing_end,
            prf=prf,
            starting_range=starting_range,
            orbit=orbit,
            beta0_coeff=beta0_coeff,
        )
        return capella_sicd


def prep_capella(granule_path: Path, work_dir: Optional[Path] = None) -> Path:
    """Prepare data for burst-based processing.

    Args:
        granule_path: Path to the UMBRA SICD file
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    reader = SICDReader(str(granule_path))
    sicd = reader.get_sicds_as_tuple()[0]
    capella_sicd = CapellaSICD.from_sarpy_sicd(sicd, granule_path)
    dem_path = work_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, capella_sicd.footprint)
    return capella_sicd, dem_path
