import argparse
from pathlib import Path
from typing import Optional

import s1reader

from multirtc.create_rtc import run_single_job
from multirtc.define_geogrid import generate_geogrids
from multirtc.prep_burst import prep_burst
from multirtc.rtc_options import RtcOptions


def opera_rtc_s1_burst(
    granule: str,
    resolution: int = 30,
    work_dir: Optional[Path] = None,
) -> Path:
    """Prepare data for SLC-based processing.

    Args:
        granules: List of Sentinel-1 level-1 granules to back-project
        resolution: Resolution of the output RTC (m)
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    input_dir = work_dir / 'input'
    output_dir = work_dir / 'output'
    [d.mkdir(parents=True, exist_ok=True) for d in [input_dir, output_dir]]

    granule_path, orbit_path, db_path, dem_path = prep_burst([granule], work_dir=input_dir)
    burst = s1reader.load_bursts(str(granule_path), str(orbit_path), 1, 'VV')[0]
    opts = RtcOptions(
        dem_path=str(dem_path),
        output_dir=str(output_dir),
        x_spacing=resolution,
        y_spacing=resolution,
    )
    geogrid = generate_geogrids(burst, opts.x_spacing, opts.y_spacing, opts.x_snap, opts.y_snap)
    run_single_job(granule, burst, geogrid, opts)


def main():
    """Create an OPERA RTC for an SLC granule

    Example command:
    multirtc S1_245714_IW1_20240809T141633_VV_6B31-BURST --resolution 40
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('granule', help='Data granule to create an RTC for.')
    parser.add_argument('--resolution', default=30, type=int, help='Resolution of the output RTC (m)')
    parser.add_argument('--work-dir', type=Path, default=None, help='Working directory for processing')
    args = parser.parse_args()

    if args.granule.endswith('-BURST'):
        opera_rtc_s1_burst(**args.__dict__)
    else:
        raise NotImplementedError('Only Sentinel-1 burst processing is supported at this time')


if __name__ == '__main__':
    main()
