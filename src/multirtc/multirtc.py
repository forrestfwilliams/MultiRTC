import argparse
from pathlib import Path
from typing import Optional

import isce3
from burst2safe.burst2safe import burst2safe

from multirtc import dem, orbit
from multirtc.create_rtc import pfa_prototype_geocode, rtc
from multirtc.prep_burst import S1BurstSlc
from multirtc.rtc_options import RtcOptions
from multirtc.sicd import SicdPfaSlc, SicdRzdSlc


def print_wkt(sicd):
    radar_grid = sicd.as_isce3_radargrid()
    dem = isce3.geometry.DEMInterpolator(sicd.hae)
    doppler = sicd.get_doppler_centroid_grid()
    wkt = isce3.geometry.get_geo_perimeter_wkt(
        grid=radar_grid, orbit=sicd.orbit, doppler=doppler, dem=dem, points_per_edge=3
    )
    print(wkt)


def prep_dirs(work_dir: Optional[Path] = None) -> tuple[Path, Path]:
    if work_dir is None:
        work_dir = Path.cwd()
    input_dir = work_dir / 'input'
    output_dir = work_dir / 'output'
    [d.mkdir(parents=True, exist_ok=True) for d in [input_dir, output_dir]]
    return input_dir, output_dir


def run_multirtc(platform: str, granule: str, resolution: int, work_dir: Path) -> None:
    """Create an OPERA RTC"""
    input_dir, output_dir = prep_dirs(work_dir)
    if platform == 'S1':
        safe_path = burst2safe(granules=[granule], all_anns=True, work_dir=input_dir)
        orbit_path = orbit.get_orbit(safe_path.with_suffix('').name, save_dir=input_dir)
        slc = S1BurstSlc(safe_path, orbit_path, granule)
    elif platform in ['CAPELLA', 'UMBRA']:
        sicd_class = {'CAPELLA': SicdRzdSlc, 'UMBRA': SicdPfaSlc}[platform]
        granule_path = input_dir / granule
        if not granule_path.exists():
            raise FileNotFoundError(f'SICD must be present in input dir {input_dir} for processing.')
        slc = sicd_class(granule_path)
    else:
        raise ValueError(f'Unsupported platform {platform}. Supported platforms are S1, CAPELLA, UMBRA.')

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
    """Create an OPERA RTC for a multiple satellite platforms

    Example command:
    multirtc UMBRA umbra_image.ntif --resolution 40
    """
    supported = ['S1', 'UMBRA', 'CAPELLA']
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('platform', choices=supported, help='Platform to create RTC for')
    parser.add_argument('granule', help='Data granule to create an RTC for.')
    parser.add_argument('--resolution', default=30, type=float, help='Resolution of the output RTC (m)')
    parser.add_argument('--work-dir', type=Path, default=None, help='Working directory for processing')
    args = parser.parse_args()

    if args.work_dir is None:
        args.work_dir = Path.cwd()
    run_multirtc(args.platform, args.granule, args.resolution, args.work_dir)


if __name__ == '__main__':
    main()
