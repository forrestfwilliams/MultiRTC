"""Create a geocoded dataset for a multiple satellite platforms"""

from pathlib import Path

from multirtc.multirtc import DOWNLOAD_SOURCES, SUPPORTED, run_multirtc


def create_parser(parser):
    parser.add_argument('platform', choices=SUPPORTED, help='Platform to create geocoded dataset for')
    parser.add_argument('granule', help='Data granule to create geocoded dataset for.')
    parser.add_argument('--resolution', type=float, help='Resolution of the output dataset (m)')
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
        apply_rtc=False,
        download_source=args.download_source,
    )
