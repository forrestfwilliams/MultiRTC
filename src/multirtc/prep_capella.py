from pathlib import Path
from typing import Optional

from osgeo import gdal

from multirtc import dem
from multirtc.sicd import SicdRzdSlc


gdal.UseExceptions()


def prep_capella(granule_path: Path, work_dir: Optional[Path] = None) -> Path:
    """Prepare data for burst-based processing.

    Args:
        granule_path: Path to the UMBRA SICD file
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    capella_sicd = SicdRzdSlc(granule_path)
    dem_path = work_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, capella_sicd.footprint)
    return capella_sicd, dem_path
