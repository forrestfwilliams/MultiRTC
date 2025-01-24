from pathlib import Path
from typing import Union
from zipfile import ZipFile

import earthaccess
import lxml.etree as ET
import requests
from shapely.geometry import Polygon, box


def download_file(
    url: str,
    download_path: Union[Path, str] = '.',
    chunk_size=10 * (2**20),
) -> Path:
    """Download a file without authentication.

    Args:
        url: URL of the file to download
        download_path: Path to save the downloaded file to
        chunk_size: Size to chunk the download into

    Returns:
        download_path: The path to the downloaded file
    """
    session = requests.Session()

    with session.get(url, stream=True) as s:
        s.raise_for_status()
        with open(download_path, 'wb') as f:
            for chunk in s.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    session.close()


def download_burst_db(save_dir: Path) -> Path:
    """Download the OPERA burst database.
    Currently using a version created using opera-adt/burst_db v0.4.0, but hope to switch to ASF-provide source.

    Args:
        save_dir: Directory to save the database to

    Returns:
        Path to the downloaded database
    """
    db_path = save_dir / 'opera-burst-bbox-only.sqlite3'

    if db_path.exists():
        return db_path

    url = 'https://ffwilliams2-shenanigans.s3.us-west-2.amazonaws.com/opera/opera-burst-bbox-only.sqlite3'
    download_file(url, db_path)
    return db_path


def download_s1_granule(granule, save_dir: Path) -> Path:
    out_path = save_dir / f'{granule}.zip'
    if out_path.exists():
        return out_path

    mission = granule[0] + granule[2]
    product_type = granule[7:10]
    url = f'https://sentinel1.asf.alaska.edu/{product_type}/{mission}/{granule}.zip'
    earthaccess.download(url, str(save_dir))
    return out_path


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
