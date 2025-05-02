from pathlib import Path
from shutil import make_archive
from typing import Optional
from zipfile import ZipFile

import lxml.etree as ET
import s1reader
from burst2safe.burst2safe import burst2safe
from shapely.geometry import Polygon, box

from multirtc import dem, orbit


def get_s1_granule_bbox(granule_path: Path, buffer: float = 0.025) -> box:
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
    footprint = Polygon(coords).buffer(buffer)
    return box(*footprint.bounds)


def prep_burst(granule: str, work_dir: Optional[Path] = None) -> Path:
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
        granule_path = burst2safe(granules=[granule], all_anns=True, work_dir=work_dir)
        make_archive(base_name=str(granule_path.with_suffix('')), format='zip', base_dir=str(granule_path))
        granule = granule_path.with_suffix('').name
        granule_path = granule_path.with_suffix('.zip')
    else:
        granule_path = work_dir / list(work_dir.glob('S1*.zip'))[0].name

    orbit_path = orbit.get_orbit(granule_path.with_suffix('').name, save_dir=work_dir)

    dem_path = work_dir / 'dem.tif'
    granule_bbox = get_s1_granule_bbox(granule_path)
    dem.download_opera_dem_for_footprint(dem_path, granule_bbox)
    burst = s1reader.load_bursts(str(granule_path), str(orbit_path), 1, 'VV')[0]
    return burst, dem_path
