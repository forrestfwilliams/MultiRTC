from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from osgeo import gdal
from skimage.registration import phase_cross_correlation


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


def get_tiling_schema(reference_path, secondary_path, tile_size=512):
    ref_trans, ref_bounds = get_geo_info(reference_path)
    sec_trans, sec_bounds = get_geo_info(secondary_path)
    inter_bounds = (
        max(ref_bounds[0], sec_bounds[0]),
        min(ref_bounds[1], sec_bounds[1]),
        max(ref_bounds[2], sec_bounds[2]),
        min(ref_bounds[3], sec_bounds[3]),
    )
    row_offset = int((inter_bounds[3] - ref_bounds[3]) / np.abs(ref_trans[5]))
    col_offset = int((inter_bounds[0] - ref_trans[0]) / np.abs(ref_trans[1]))
    nrows = (inter_bounds[3] - inter_bounds[1]) / np.abs(ref_trans[5])
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


def rle(reference_path: Path, secondary_path: Path, project: str, basedir: Path):
    project_dir = basedir / project
    project_dir.mkdir(parents=True, exist_ok=True)
    min_val, max_val = get_flattened_range(reference_path)
    pixel_size = get_geo_info(reference_path)[0][1]
    tiles = get_tiling_schema(reference_path, secondary_path)
    rows = []
    for tile in tiles:
        ref_data, sec_data = load_tiles(reference_path, secondary_path, tile, val_bounnds=(min_val, max_val))
        n_pixels = np.prod(tile.shape)
        if np.isnan(ref_data).sum() > 0.1 * n_pixels or np.isnan(sec_data).sum() > 0.1 * n_pixels:
            continue
        ref_data[np.isnan(ref_data)] = 0.0
        sec_data[np.isnan(sec_data)] = 0.0
        shift, phase, error = phase_cross_correlation(ref_data, sec_data, upsample_factor=10)
        row = {
            'id': tile.id,
            'shift_x': shift[1] * pixel_size,
            'shift_y': shift[0] * pixel_size,
            'error': error * pixel_size,
            'phase': phase,
        }
        rows.append(pd.Series(row))
    df = pd.DataFrame(rows)
    df.to_csv(project_dir / 'rle_results.csv', index=False)


def create_parser(parser):
    parser.add_argument('reference', type=str, help='Path to the reference image')
    parser.add_argument('secondary', type=str, help='Path to the secondary image')
    parser.add_argument('project', type=str, help='Directory to save the results')
    parser.add_argument('--basedir', type=str, default='.', help='Base directory for the project')
    return parser


def run(args):
    args.reference = Path(args.reference)
    assert args.reference.exists(), f'RTC file {args.reference} does not exist'
    args.secondary = Path(args.secondary)
    assert args.secondary.exists(), f'RTC file {args.secondary} does not exist'
    args.basedir = Path(args.basedir)
    assert args.basedir.exists(), f'Base directory {args.basedir} does not exist'
    rle(args.reference, args.secondary, project=args.project, basedir=args.basedir)


def main():
    parser = ArgumentParser(description='Relative Location Error (RLE) analysis')
    parser = create_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
