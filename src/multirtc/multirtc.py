import argparse
from pathlib import Path
from typing import Optional

from burst2safe.burst2safe import burst2safe
from s1reader.s1_orbit import retrieve_orbit_file

from multirtc import dem
from multirtc.base import Slc
from multirtc.create_rtc import pfa_prototype_geocode, rtc
from multirtc.rtc_options import RtcOptions
from multirtc.sentinel1 import S1BurstSlc
from multirtc.sicd import SicdPfaSlc, SicdRzdSlc


SUPPORTED = ['S1', 'UMBRA', 'CAPELLA', 'ICEYE']


def prep_dirs(work_dir: Optional[Path] = None) -> tuple[Path, Path]:
    """Prepare input and output directories for processing.

    Args:
        work_dir: Working directory. If None, current working directory is used.

    Returns:
        Tuple of input and output directories.
    """
    if work_dir is None:
        work_dir = Path.cwd()
    input_dir = work_dir / 'input'
    output_dir = work_dir / 'output'
    [d.mkdir(parents=True, exist_ok=True) for d in [input_dir, output_dir]]
    return input_dir, output_dir


def get_slc(platform: str, granule: str, input_dir: Path) -> Slc:
    """
    Get the SLC object for the specified platform and granule.

    Args:
        platform: Platform type (e.g., 'UMBRA').
        granule: Granule name if data is available in ASF archive, or filename if granule is already downloaded.
        input_dir: Directory containing the input data.

    Returns:
        Slc subclass object for the specified platform and granule.
    """
    if platform == 'S1':
        safe_path = burst2safe(granules=[granule], all_anns=True, work_dir=input_dir)
        orbit_path = Path(retrieve_orbit_file(safe_path.name, str(input_dir), concatenate=True))
        slc = S1BurstSlc(safe_path, orbit_path, granule)
    elif platform in ['CAPELLA', 'ICEYE', 'UMBRA']:
        sicd_class = {'CAPELLA': SicdRzdSlc, 'ICEYE': SicdRzdSlc, 'UMBRA': SicdPfaSlc}[platform]
        granule_path = input_dir / granule
        if not granule_path.exists():
            raise FileNotFoundError(f'SICD must be present in input dir {input_dir} for processing.')
        slc = sicd_class(granule_path)
    else:
        raise ValueError(f'Unsupported platform {platform}. Supported platforms are {",".join(SUPPORTED)}.')
    return slc


def run_multirtc(platform: str, granule: str, resolution: int, work_dir: Path) -> None:
    """Create an RTC or Geocoded dataset using the OPERA algorithm.

    Args:
        platform: Platform type (e.g., 'UMBRA').
        granule: Granule name if data is available in ASF archive, or filename if granule is already downloaded.
        resolution: Resolution of the output RTC (in meters).
        work_dir: Working directory for processing.
    """
    input_dir, output_dir = prep_dirs(work_dir)
    slc = get_slc(platform, granule, input_dir)
    dem_path = input_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, slc.footprint)
    geogrid = slc.create_geogrid(spacing_meters=resolution)
    if slc.supports_rtc:
        opts = RtcOptions(
            dem_path=str(dem_path),
            output_dir=str(output_dir),
            resolution=resolution,
            apply_bistatic_delay=slc.supports_bistatic_delay,
            apply_static_tropo=slc.supports_static_tropo,
        )
        rtc(slc, geogrid, opts)
    else:
        pfa_prototype_geocode(slc, geogrid, dem_path, output_dir)


def main():
    """Create a RTC or geocoded dataset for a multiple satellite platforms

    Example command:
    multirtc UMBRA umbra_image.ntif --resolution 40
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('platform', choices=SUPPORTED, help='Platform to create RTC for')
    parser.add_argument('granule', help='Data granule to create an RTC for.')
    parser.add_argument('--resolution', default=30, type=float, help='Resolution of the output RTC (m)')
    parser.add_argument('--work-dir', type=Path, default=None, help='Working directory for processing')
    args = parser.parse_args()

    if args.work_dir is None:
        args.work_dir = Path.cwd()
    run_multirtc(args.platform, args.granule, args.resolution, args.work_dir)


if __name__ == '__main__':
    main()
