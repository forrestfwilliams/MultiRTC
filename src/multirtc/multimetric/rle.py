"""Relative Location Error (RLE) analysis"""

from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from osgeo import gdal
from scipy.stats import median_abs_deviation
from skimage.registration import phase_cross_correlation
from tqdm import tqdm


gdal.UseExceptions()


@dataclass
class Tile:
    id: str
    ref_row: int
    ref_col: int
    shape: tuple
    bounds: tuple


def get_flattened_range(image_path, pct_min=1, pct_max=99):
    ds = gdal.Open(str(image_path), gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    return np.nanpercentile(data, pct_min), np.nanpercentile(data, pct_max)


def get_geo_info(image_path):
    ds = gdal.Open(str(image_path), gdal.GA_ReadOnly)
    trans = ds.GetGeoTransform()
    nrow, ncol = ds.RasterYSize, ds.RasterXSize
    ds = None
    bounds = (
        trans[0],
        trans[3] + trans[5] * nrow,
        trans[0] + trans[1] * ncol,
        trans[3],
    )
    return trans, bounds


def get_tiling_schema(reference_path, secondary_path, tile_size=1024):
    ref_trans, ref_bounds = get_geo_info(reference_path)
    sec_trans, sec_bounds = get_geo_info(secondary_path)
    inter_bounds = (
        max(ref_bounds[0], sec_bounds[0]),
        max(ref_bounds[1], sec_bounds[1]),
        min(ref_bounds[2], sec_bounds[2]),
        min(ref_bounds[3], sec_bounds[3]),
    )
    row_offset = max(int((inter_bounds[3] - ref_bounds[3]) / np.abs(ref_trans[5])), tile_size)
    col_offset = max(int((inter_bounds[0] - ref_bounds[0]) / np.abs(ref_trans[1])), 0)
    nrows = ((inter_bounds[3] - inter_bounds[1]) / np.abs(ref_trans[5])) - tile_size
    nrow_tiles = int(np.floor(nrows / tile_size))
    ncols = (inter_bounds[2] - inter_bounds[0]) / np.abs(ref_trans[1])
    ncol_tiles = int(np.floor(ncols / tile_size))
    tiles = []
    for irow, icol in np.ndindex(nrow_tiles, ncol_tiles):
        minrow = row_offset + (irow * tile_size)
        mincol = col_offset + (icol * tile_size)
        maxy = ref_trans[3] + (minrow * ref_trans[5])
        miny = maxy + (tile_size * ref_trans[5])
        minx = ref_trans[0] + (mincol * ref_trans[1])
        maxx = minx + (tile_size * ref_trans[1])
        bounds = (minx, miny, maxx, maxy)
        tile = Tile(f'tile_{irow}_{icol}', minrow, mincol, (tile_size, tile_size), bounds)
        tiles.append(tile)
    return tiles


def load_tiles(reference_path: Path, secondary_path: Path, tile: Tile, val_bounnds=None):
    ref_ds = gdal.Open(str(reference_path), gdal.GA_ReadOnly)
    ref_band = ref_ds.GetRasterBand(1)
    ref_data = ref_band.ReadAsArray(tile.ref_col, tile.ref_row, *tile.shape)
    with NamedTemporaryFile(suffix='.tif') as sec_warped:
        gdal.Warp(
            sec_warped.name,
            str(secondary_path),
            outputBounds=tile.bounds,
            width=tile.shape[1],
            height=tile.shape[0],
            resampleAlg=gdal.GRA_Bilinear,
            format='GTiff',
        )
        sec_ds = gdal.Open(sec_warped.name, gdal.GA_ReadOnly)
        sec_band = sec_ds.GetRasterBand(1)
        sec_data = sec_band.ReadAsArray()
        sec_ds = None
        assert sec_data.shape == ref_data.shape, 'Reference and secondary tile shapes do not match'

    if val_bounnds is not None:
        ref_data[ref_data <= val_bounnds[0]] = np.nan
        ref_data[ref_data >= val_bounnds[1]] = np.nan
        sec_data[sec_data <= val_bounnds[0]] = np.nan
        sec_data[sec_data >= val_bounnds[1]] = np.nan
    return ref_data, sec_data


def plot_offsets(df: pd.DataFrame, output_path: Path):
    row_col = [id.split('_')[1:] for id in df['id']]
    rows = [int(row) for row, _ in row_col]
    minrow, maxrow = min(rows), max(rows)
    cols = [int(col) for _, col in row_col]
    mincol, maxcol = min(cols), max(cols)
    blank = np.zeros((maxrow - minrow + 1, maxcol - mincol + 1))
    blank[:, :] = np.nan
    y_offset = blank.copy()
    x_offset = blank.copy()
    for idx, row in df.iterrows():
        if not row['valid']:
            continue
        i, j = map(int, row['id'].split('_')[1:])
        y_offset[i - minrow, j - mincol] = row['shift_y']
        x_offset[i - minrow, j - mincol] = row['shift_x']
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    vmin = np.nanmean(y_offset.flatten()) - 3 * np.nanstd(y_offset.flatten())
    vmax = np.nanmean(y_offset.flatten()) + 3 * np.nanstd(y_offset.flatten())
    im1 = ax1.imshow(y_offset, cmap='coolwarm', vmin=vmin, vmax=vmax)
    ax1.xaxis.set_ticks_position('top')  # Moves the ticks
    ax1.xaxis.set_label_position('top')
    ax1.tick_params(axis='x', bottom=False)
    ax1.set_title('Y Offset (m)')
    ax1.set_xlabel('Column Index')
    ax1.set_ylabel('Row Index')
    vmin = np.nanmean(x_offset.flatten()) - 3 * np.nanstd(x_offset.flatten())
    vmax = np.nanmean(x_offset.flatten()) + 3 * np.nanstd(x_offset.flatten())
    im2 = ax2.imshow(x_offset, cmap='coolwarm', vmin=vmin, vmax=vmax)
    ax2.xaxis.set_ticks_position('top')  # Moves the ticks
    ax2.xaxis.set_label_position('top')
    ax2.tick_params(axis='x', bottom=False)
    ax2.set_title('X Offset (m)')
    ax2.set_xlabel('Column Index')
    ax2.set_ylabel('Row Index')
    plt.colorbar(im1, ax=ax1)
    plt.colorbar(im2, ax=ax2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')


def rle(reference_path: Path, secondary_path: Path, project: str, basedir: Path, max_nan_ratio=0.1, max_error=0.01):
    project_dir = basedir / project
    project_dir.mkdir(parents=True, exist_ok=True)
    min_val, max_val = get_flattened_range(reference_path)
    pixel_size = get_geo_info(reference_path)[0][1]
    tiles = get_tiling_schema(reference_path, secondary_path)
    rows = []
    for tile in tqdm(tiles, desc='Processing tiles'):
        ref_data, sec_data = load_tiles(reference_path, secondary_path, tile, val_bounnds=(min_val, max_val))
        max_nans = max_nan_ratio * np.prod(tile.shape)
        if np.isnan(ref_data).sum() > max_nans or np.isnan(sec_data).sum() > max_nans:
            continue
        ref_data[np.isnan(ref_data)] = 0.0
        sec_data[np.isnan(sec_data)] = 0.0
        shift, _, error = phase_cross_correlation(ref_data, sec_data, upsample_factor=16, overlap_ratio=0.8)
        shift_y, shift_x = shift * pixel_size
        rows.append(pd.Series({'id': tile.id, 'shift_x': shift_x, 'shift_y': shift_y, 'error': error}))
    base_name = f'{reference_path.stem}_x_{secondary_path.stem}'
    df = pd.DataFrame(rows)

    df['shift_x_mdz'] = 0.6745 * np.abs(df['shift_x'] - np.median(df['shift_x'])) / median_abs_deviation(df['shift_x'])
    df['shift_y_mdz'] = 0.6745 * np.abs(df['shift_y'] - np.median(df['shift_y'])) / median_abs_deviation(df['shift_y'])
    df['valid'] = (df['shift_x_mdz'] <= 3.5) & (df['shift_y_mdz'] <= 3.5)
    plot_offsets(df, project_dir / f'{base_name}_offsets.png')
    mean_row = df.loc[df['valid'], df.columns.drop(['id', 'valid'])].mean()
    mean_row['id'] = 'mean'
    std_row = df.loc[df['valid'], df.columns.drop(['id', 'valid'])].std()
    std_row['id'] = 'std'
    df = pd.concat([df, pd.DataFrame([mean_row, std_row])], ignore_index=True)
    df.to_csv(project_dir / f'{base_name}.csv', index=False)


def create_parser(parser):
    parser.add_argument('reference', type=str, help='Path to the reference image')
    parser.add_argument('secondary', type=str, help='Path to the secondary image')
    parser.add_argument('project', type=str, help='Directory to save the results')
    parser.add_argument('--basedir', type=str, default='.', help='Base directory for the project')
    return parser


def run(args):
    args.reference = Path(args.reference)
    assert args.reference.exists(), f'Image file {args.reference} does not exist'
    args.secondary = Path(args.secondary)
    assert args.secondary.exists(), f'Image file {args.secondary} does not exist'
    args.basedir = Path(args.basedir)
    assert args.basedir.exists(), f'Base directory {args.basedir} does not exist'
    rle(args.reference, args.secondary, project=args.project, basedir=args.basedir)


def main():
    parser = ArgumentParser(description=__doc__)
    parser = create_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
