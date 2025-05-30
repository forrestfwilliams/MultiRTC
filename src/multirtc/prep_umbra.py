from pathlib import Path
from typing import Optional

from multirtc import dem
from multirtc.sicd import SicdPfaSlc


def prep_umbra(granule_path: Path, work_dir: Optional[Path] = None) -> Path:
    """Prepare data for burst-based processing.

    Args:
        granule_path: Path to the UMBRA SICD file
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    umbra_sicd = SicdPfaSlc(granule_path)
    dem_path = work_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, umbra_sicd.footprint)
    return umbra_sicd, dem_path
