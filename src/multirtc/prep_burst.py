from pathlib import Path
from shutil import make_archive
from typing import Optional
from zipfile import ZipFile

import isce3
import lxml.etree as ET
import s1reader
from burst2safe.burst2safe import burst2safe
from shapely.geometry import Polygon, box

from multirtc import dem, orbit
from multirtc.base import SlcTemplate, from_isce_datetime, to_isce_datetime


class S1BurstSlc(SlcTemplate):
    def __init__(self, safe_path, orbit_path, burst_name):
        _, burst_id, swath, _, polarization, _ = burst_name.split('_')
        burst_id = int(burst_id)
        swath_num = int(swath[2])
        bursts = s1reader.load_bursts(str(safe_path), str(orbit_path), swath_num, polarization)
        burst = [b for b in bursts if str(b.burst_id).endswith(f'{burst_id}_{swath.lower()}')][0]
        del bursts
        vrt_path = safe_path.parent / f'{burst_name}.vrt'
        burst.slc_to_vrt_file(vrt_path)
        self.id = burst_name
        self.filepath = vrt_path
        self.footprint = burst.border[0]
        self.center = burst.center
        self.lookside = 'right'
        self.wavelength = burst.wavelength
        self.polarization = burst.polarization
        self.shape = burst.shape
        self.range_pixel_spacing = burst.range_pixel_spacing
        self.reference_time = from_isce_datetime(burst.orbit.reference_epoch)
        self.sensing_start = (burst.sensing_start - self.reference_time).total_seconds()
        self.starting_range = burst.starting_range
        self.prf = 1 / burst.azimuth_time_interval
        self.orbit = burst.orbit
        self.doppler_centroid_grid = isce3.core.LUT2d()
        self.radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start,
            wavelength=self.wavelength,
            prf=self.prf,
            starting_range=self.starting_range,
            range_pixel_spacing=self.range_pixel_spacing,
            lookside=isce3.core.LookSide.Right,
            length=self.shape[0],
            width=self.shape[1],
            ref_epoch=to_isce_datetime(self.reference_time),
        )
        self.first_valid_line = burst.first_valid_line
        self.last_valid_line = burst.last_valid_line
        self.first_valid_sample = burst.first_valid_sample
        self.last_valid_sample = burst.last_valid_sample
        self.source = burst


def get_s1_granule_bbox(granule_path: Path) -> box:
    if granule_path.suffix == '.zip':
        with ZipFile(granule_path, 'r') as z:
            manifest_path = [x for x in z.namelist() if x.endswith('manifest.safe')][0]
            with z.open(manifest_path) as m:
                manifest = ET.parse(m).getroot()
    else:
        manifest_path = granule_path / 'manifest.safe'
        manifest = ET.parse(manifest_path).getroot()

    frame_element = [x for x in manifest.findall('.//metadataObject') if x.get('ID') == 'measurementFrameSet'][0]
    frame_string = frame_element.find('.//{http://www.opengis.net/gml}coordinates').text
    coord_strings = [pair.split(',') for pair in frame_string.split(' ')]
    coords = [(float(lon), float(lat)) for lat, lon in coord_strings]
    return Polygon(coords)


def prep_burst(burst_granule: str, work_dir: Optional[Path] = None) -> Path:
    """Prepare data for burst-based processing.

    Args:
        granule: Sentinel-1 burst SLC granule to create RTC dataset for
        use_resorb: Use the RESORB orbits instead of the POEORB orbits
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()

    print('Downloading data...')

    if len(list(work_dir.glob('S1*.zip'))) == 0:
        granule_path = burst2safe(granules=[burst_granule], all_anns=True, work_dir=work_dir)
        make_archive(base_name=str(granule_path.with_suffix('')), format='zip', base_dir=str(granule_path))
        granule_path = granule_path.with_suffix('.zip')
    else:
        granule_path = work_dir / list(work_dir.glob('S1*.zip'))[0].name

    orbit_path = orbit.get_orbit(granule_path.with_suffix('').name, save_dir=work_dir)

    burst_slc = S1BurstSlc(granule_path, orbit_path, burst_granule)
    dem_path = work_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, burst_slc.footprint)
    return burst_slc, dem_path
