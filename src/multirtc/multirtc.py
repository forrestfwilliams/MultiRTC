"""Create an RTC dataset for a multiple satellite platforms"""

import logging
from pathlib import Path

from burst2safe.burst2safe import burst2safe
from s1reader.s1_orbit import retrieve_orbit_file

from multirtc import dem
from multirtc.base import Slc
from multirtc.cdse import burst_to_parent_slc, download_slc_from_cdse, ensure_cdse_credentials
from multirtc.create_rtc import rtc
from multirtc.rtc_options import RtcOptions
from multirtc.sentinel1 import S1BurstSlc
from multirtc.sicd import SicdPfaSlc, SicdRzdSlc


logger = logging.getLogger(__name__)

DOWNLOAD_SOURCES = ['ASF', 'CDSE']


SUPPORTED = ['S1', 'UMBRA', 'CAPELLA', 'ICEYE']


def prep_dirs(work_dir: Path | None = None) -> tuple[Path, Path]:
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


def get_slc(platform: str, granule: str, input_dir: Path, download_source: str = 'ASF') -> Slc:
    """
    Get the SLC object for the specified platform and granule.

    Args:
        platform: Platform type (e.g., 'UMBRA').
        granule: Granule name if data is available in ASF archive, or filename if granule is already downloaded.
        input_dir: Directory containing the input data.
        download_source: Source for downloading Sentinel-1 SLC data ('ASF' or 'CDSE').

    Returns:
        Slc subclass object for the specified platform and granule.
    """
    if platform == 'S1':
        if download_source == 'CDSE':
            ensure_cdse_credentials()
            parent_slc = burst_to_parent_slc(granule)
            safe_path = download_slc_from_cdse(parent_slc, input_dir)
        else:
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


def run_multirtc(
    platform: str,
    granule: str,
    resolution: int,
    work_dir: Path,
    dem_path: Path | None = None,
    apply_rtc=True,
    download_source: str = 'ASF',
) -> None:
    """Create an RTC or Geocoded dataset using the OPERA algorithm.

    Args:
        platform: Platform type (e.g., 'UMBRA').
        granule: Granule name if data is available in ASF archive, or filename if granule is already downloaded.
        resolution: Resolution of the output RTC (in meters).
        work_dir: Working directory for processing.
        dem_path: Path to the DEM to use for processing. If None, the NISAR DEM will be downloaded.
        apply_rtc: If True perform radiometric correction; if False, only geocode.
        download_source: Source for downloading Sentinel-1 SLC data ('ASF' or 'CDSE').
    """
    input_dir, output_dir = prep_dirs(work_dir)
    slc = get_slc(platform, granule, input_dir, download_source=download_source)
    if dem_path is None:
        dem_path = input_dir / 'dem.tif'
        dem.download_opera_dem_for_footprint(dem_path, slc.footprint)
    dem.validate_dem(dem_path, slc.footprint)
    geogrid = slc.create_geogrid(spacing_meters=resolution)
    if slc.supports_rtc:
        opts = RtcOptions(
            dem_path=str(dem_path),
            output_dir=str(output_dir),
            apply_rtc=apply_rtc,
            resolution=resolution,
            apply_bistatic_delay=slc.supports_bistatic_delay,
            apply_static_tropo=slc.supports_static_tropo,
        )
        rtc(slc, geogrid, opts)
    else:
        raise NotImplementedError(
            'RTC creation is not supported for this input. For polar grid support, use the multirtc docker image:\n'
            'https://github.com/forrestfwilliams/MultiRTC/pkgs/container/multirtc'
        )


def create_parser(parser):
    parser.add_argument('platform', choices=SUPPORTED, help='Platform to create RTC for')
    parser.add_argument('granule', help='Data granule to create an RTC for.')
    parser.add_argument('--resolution', type=float, help='Resolution of the output RTC (m)')
    parser.add_argument('--dem', type=Path, default=None, help='Path to the DEM to use for processing')
    parser.add_argument('--work-dir', type=Path, default=None, help='Working directory for processing')
    parser.add_argument(
        '--download-source',
        type=str,
        choices=DOWNLOAD_SOURCES,
        default='ASF',
        help=(
            "Source for downloading Sentinel-1 SLC data. "
            "'ASF' uses Alaska Satellite Facility (default). "
            "'CDSE' uses Copernicus Data Space Ecosystem "
            "(requires CDSE_USERNAME/CDSE_PASSWORD env vars or ~/.netrc)."
        ),
    )
    return parser


def run(args):
    if args.dem is not None:
        assert args.dem.exists(), f'DEM file {args.dem} does not exist.'
    if args.work_dir is None:
        args.work_dir = Path.cwd()
    run_multirtc(
        args.platform,
        args.granule,
        args.resolution,
        args.work_dir,
        args.dem,
        apply_rtc=True,
        download_source=args.download_source,
    )
