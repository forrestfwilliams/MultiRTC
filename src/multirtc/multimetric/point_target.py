"""Point target (resolution, PSLR, and ISLR) analysis"""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from isce3.cal.point_target_info import analyze_point_target_chip, dB
from shapely.geometry import box

from multirtc.multimetric import corner_reflector
from multirtc.multirtc import get_slc
from multirtc.sicd import reformat_for_isce


SUPPORTED = ['UMBRA']


def plot_profile(ax, n, magnitude, phase, title=None):
    peak = abs(magnitude).max()
    ax.plot(n, dB(magnitude) - dB(peak), '-k')
    ax.set_ylim((-40, 0.3))
    ax.set_ylabel('Power (dB)')
    ax_twin = ax.twinx()
    phase_color = '0.75'
    ax_twin.plot(n, phase, color=phase_color)
    ax_twin.set_ylim((-np.pi, np.pi))
    ax_twin.set_ylabel('Phase (rad)')
    ax_twin.spines['right'].set_color(phase_color)
    ax_twin.tick_params(axis='y', colors=phase_color)
    ax_twin.yaxis.label.set_color(phase_color)
    ax.set_xlim((-15, 15))
    ax.spines['top'].set_visible(False)
    ax_twin.spines['top'].set_visible(False)
    ax.set_title(title)
    return ax


def analyze_point_targets(platform, filepath, project, basedir, width=64):
    outdir = Path(basedir).expanduser() / project
    outdir.mkdir(exist_ok=True, parents=True)
    slc = get_slc('UMBRA', filepath.name, filepath.parent.expanduser())
    cr_df = corner_reflector.get_cr_df(box(*slc.footprint.bounds), 4326, slc.reference_time, slc.look_angle, outdir)
    cr_df = corner_reflector.add_rdr_image_location(slc, cr_df, search_radius=50)
    blank = [np.nan] * cr_df.shape[0]
    cr_df = cr_df.assign(az_res=blank, az_pslr=blank, az_islr=blank, rng_res=blank, rng_pslr=blank, rng_islr=blank)
    for idx, row in cr_df.iterrows():
        center = (int(row['yloc']), int(row['xloc']))
        half_width = int(np.floor(width / 2))
        row_start = center[0] - half_width
        row_end = center[0] + half_width
        col_start = center[1] - half_width
        col_end = center[1] + half_width
        chip = slc.load_scaled_data('beta0', rowrange=(row_start, row_end), colrange=(col_start, col_end))
        chip = reformat_for_isce(chip, slc.az_reversed)
        data, _ = analyze_point_target_chip(
            chip, chip_min_i=row_start, chip_min_j=col_start, i_pos=center[0], j_pos=center[1], cuts=True
        )
        row['az_res'] = data['azimuth']['resolution']
        row['az_pslr'] = data['azimuth']['PSLR']
        row['az_islr'] = data['azimuth']['ISLR']
        row['rng_res'] = data['range']['resolution']
        row['rng_pslr'] = data['range']['PSLR']
        row['rng_islr'] = data['range']['ISLR']
        cr_df.iloc[idx] = row

        figure, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
        peak = abs(chip).max()
        ax1.imshow(dB(chip) - dB(peak), cmap='gray', vmin=-40, vmax=0, interpolation='none')
        ax1.set_title(f'CR {int(row["ID"])}')
        ax1.set_xlabel('Range')
        ax1.set_ylabel('Azimuth')
        ax1.axis('off')
        ax2 = plot_profile(
            ax2,
            np.array(data['azimuth']['cut']),
            np.array(data['azimuth']['magnitude cut']),
            np.array(data['azimuth']['phase cut']),
            title='Azimuth Profile',
        )
        ax3 = plot_profile(
            ax3,
            np.array(data['range']['cut']),
            np.array(data['range']['magnitude cut']),
            np.array(data['range']['phase cut']),
            title='Range Profile',
        )
        plt.tight_layout()
        name_base = f'{project}_CR_{int(row["ID"])}'
        figure.savefig(outdir / f'{name_base}_point_target.png', dpi=300)

    cr_df['az_res'] = cr_df['az_res'] * slc.source.Grid.Row.SS
    cr_df['rng_res'] = cr_df['rng_res'] * slc.source.Grid.Col.SS
    cr_df.to_csv(outdir / f'{project}_point_target.csv', index=False)


def create_parser(parser):
    parser.add_argument('platform', choices=SUPPORTED, help='Platform to create RTC for')
    parser.add_argument('filepath', type=str, help='Path to the file to be processed')
    parser.add_argument('project', type=str, help='File prefix and output directory')
    parser.add_argument('--basedir', type=str, default='.', help='Base directory for the project dir')
    return parser


def run(args):
    if args.platform not in SUPPORTED:
        raise ValueError(f'Platform {args.platform} is not supported. Supported platforms: {SUPPORTED}')

    args.filepath = Path(args.filepath).expanduser()
    assert args.filepath.exists()
    args.basedir = Path(args.basedir).expanduser()
    analyze_point_targets(args.platform, args.filepath, args.project, basedir=args.basedir)


def main():
    parser = ArgumentParser(description='Point target (resolution, PSLR, and ISLR) analysis')
    parser = create_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
