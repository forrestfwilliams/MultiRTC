import argparse
from pathlib import Path
from typing import Optional

from multirtc.create_rtc import run_single_job, umbra_rtc
from multirtc.define_geogrid import generate_geogrids
from multirtc.prep_burst import prep_burst
from multirtc.prep_umbra import prep_umbra
from multirtc.rtc_options import RtcOptions


def opera_rtc_s1_burst(granule: str, resolution: int = 30, work_dir: Optional[Path] = None) -> None:
    """Create an OPERA RTC for a Sentinel-1 burst

    Args:
        granule: Sentinel-1 level-1 granule name to create an RTC for
        resolution: Resolution of the output RTC (m)
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    input_dir = work_dir / 'input'
    output_dir = work_dir / 'output'
    [d.mkdir(parents=True, exist_ok=True) for d in [input_dir, output_dir]]

    burst, dem_path = prep_burst(granule, work_dir=input_dir)
    opts = RtcOptions(dem_path=str(dem_path), output_dir=str(output_dir), resolution=resolution)
    geogrid = generate_geogrids(burst, opts.resolution)
    run_single_job(granule, burst, geogrid, opts)


def opera_rtc_umbra_sicd(granule: str, resolution: int = 30, work_dir: Optional[Path] = None) -> None:
    """Create an OPERA RTC for an UMBRA SICD file

    Args:
        granule: Umbra SICD file name to create an RTC for
        resolution: Resolution of the output RTC (m)
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    input_dir = work_dir / 'input'
    output_dir = work_dir / 'output'
    granule_path = input_dir / granule
    if not granule_path.exists():
        raise FileNotFoundError(f'Umbra SICD must be present in input dir {input_dir} for processing.')
    [d.mkdir(parents=True, exist_ok=True) for d in [input_dir, output_dir]]
    umbra_sicd, dem_path = prep_umbra(granule_path, work_dir=input_dir)
    geogrid = generate_geogrids(umbra_sicd, resolution, rda=False)
    umbra_rtc(umbra_sicd, geogrid, dem_path, output_dir=output_dir)


def main():
    """Create an OPERA RTC for an Umbra SICD SLC granule

    Example command:
    multirtc umbra_image.ntif --resolution 40
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('platform', choices=['S1', 'UMBRA'], help='Platform to create RTC for')
    parser.add_argument('granule', help='Data granule to create an RTC for.')
    parser.add_argument('--resolution', default=30, type=float, help='Resolution of the output RTC (m)')
    parser.add_argument('--work-dir', type=Path, default=None, help='Working directory for processing')
    args = parser.parse_args()

    if args.platform == 'S1':
        opera_rtc_s1_burst(args.granule, args.resolution, args.work_dir)
    elif args.platform == 'UMBRA':
        opera_rtc_umbra_sicd(args.granule, args.resolution, args.work_dir)
    else:
        raise NotImplementedError('Only Sentinel-1 burst and Umbra processing are supported at this time')


if __name__ == '__main__':
    main()
