from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import isce3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from lmfit import Model
from osgeo import gdal, osr
from pyproj import Transformer
from shapely import box
from shapely.geometry import Point
from shapely.ops import transform


gdal.UseExceptions()


def gaussfit(x, y, A, x0, y0, sigma_x, sigma_y, theta):
    theta = np.radians(theta)
    sigx2 = sigma_x**2
    sigy2 = sigma_y**2
    a = np.cos(theta) ** 2 / (2 * sigx2) + np.sin(theta) ** 2 / (2 * sigy2)
    b = np.sin(theta) ** 2 / (2 * sigx2) + np.cos(theta) ** 2 / (2 * sigy2)
    c = np.sin(2 * theta) / (4 * sigx2) - np.sin(2 * theta) / (4 * sigy2)

    expo = -a * (x - x0) ** 2 - b * (y - y0) ** 2 - 2 * c * (x - x0) * (y - y0)
    return A * np.exp(expo)


def get_cr_df(bounds, epsg, date, outdir):
    rosamond_bounds = [-124.409591, 32.534156, -114.131211, 42.009518]
    rosamond_bounds = box(*rosamond_bounds)
    transformer = Transformer.from_crs('EPSG:4326', f'EPSG:{epsg}', always_xy=True)
    rosamond_bounds = transform(transformer.transform, rosamond_bounds)
    assert bounds.intersects(rosamond_bounds), f'Images does not intersect with Rosamond bounds {rosamond_bounds}.'
    date_str = date.strftime('%Y-%m-%d+%H\u0021%M')

    crdata = outdir / f'{date_str.split("+")[0]}_crdata.csv'
    if not crdata.exists():
        res = requests.get(
            f'https://uavsar.jpl.nasa.gov/cgi-bin/corner-reflectors.pl?date={str(date_str)}&project=rosamond_plate_location'
        )
        crdata.write_bytes(res.content)
    cr_df = pd.read_csv(crdata)
    new_cols = {
        '   "Corner ID"': 'ID',
        'Latitude (deg)': 'lat',
        'Longitude (deg)': 'lon',
        'Azimuth (deg)': 'azm',
        'Height Above Ellipsoid (m)': 'hgt',
        'Side Length (m)': 'slen',
    }
    cr_df.rename(columns=new_cols, inplace=True)
    cr_df.drop(columns=cr_df.keys()[-1], inplace=True)
    return cr_df


def add_image_location(cr_df, epsg, x_start, y_start, x_spacing, y_spacing, bounds):
    blank = [np.nan] * cr_df.shape[0]
    blank_bool = [False] * cr_df.shape[0]
    cr_df = cr_df.assign(
        UTMx=blank,
        UTMy=blank,
        xloc=blank,
        yloc=blank,
        xloc_floats=blank,
        yloc_floats=blank,
        inPoly=blank_bool,
    )
    transformer = Transformer.from_crs('EPSG:4326', f'EPSG:{epsg}', always_xy=True)
    for idx, row in cr_df.iterrows():
        row['UTMx'], row['UTMy'] = transformer.transform(row['lon'], row['lat'])
        row['xloc_floats'] = (row['UTMx'] - x_start) / x_spacing
        row['xloc'] = int(round(row['xloc_floats']))
        row['yloc_floats'] = (row['UTMy'] - y_start) / y_spacing
        row['yloc'] = int(round(row['yloc_floats']))
        row['inPoly'] = bounds.intersects(Point(row['UTMx'], row['UTMy']))
        cr_df.iloc[idx] = row

    cr_df = cr_df[cr_df['inPoly']]
    cr_df.drop('inPoly', axis=1, inplace=True)
    cr_df = cr_df.reset_index(drop=True)
    return cr_df


def filter_valid_data(cr_df, data):
    cr_df = cr_df.assign(has_data=False)
    for idx, row in cr_df.iterrows():
        xloc = int(row['xloc'])
        yloc = int(row['yloc'])
        local_data = data[yloc - 2 : yloc + 2, xloc - 3 : xloc + 3]
        row['has_data'] = bool(~np.all(np.isnan(local_data)))
        cr_df.iloc[idx] = row
    cr_df = cr_df[cr_df['has_data']]
    cr_df.drop('has_data', axis=1, inplace=True)
    cr_df = cr_df.reset_index(drop=True)
    cr_df = cr_df.loc[cr_df['slen'] > 0.8].reset_index(drop=True)  # excluding SWOT CRs (0.7 m as a side length)
    return cr_df


def filter_orientation(cr_df, azimuth_angle):
    looking_east = azimuth_angle < 180
    if looking_east:
        cr_df = cr_df[(cr_df['azm'] < 200) & (cr_df['azm'] > 20)].reset_index(drop=True)
    else:
        cr_df = cr_df[cr_df['azm'] > 340].reset_index(drop=True)
    return cr_df


def plot_crs_on_image(cr_df, data, outdir, fileprefix):
    buffer = 50
    min_x = cr_df['xloc'].min() - buffer
    max_x = cr_df['xloc'].max() + buffer
    min_y = cr_df['yloc'].min() - buffer
    max_y = cr_df['yloc'].max() + buffer

    fig, ax = plt.subplots(figsize=(15, 7))
    ax.imshow(data, cmap='gray', interpolation='bilinear', vmin=0.3, vmax=1.7, origin='upper')
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.axis('off')

    for sl in pd.unique(cr_df.slen):
        xx = cr_df.loc[cr_df['slen'] == sl]['xloc']
        yy = cr_df.loc[cr_df['slen'] == sl]['yloc']
        id = cr_df.loc[cr_df['slen'] == sl]['ID']
        color = {2.4384: 'blue', 4.8: 'red', 2.8: 'yellow'}.get(sl, 'green')
        ax.scatter(xx, yy, color=color, marker='o', facecolor='none', lw=1)
        for id_ann, xx_ann, yy_ann in zip(id, xx, yy):
            id_ann = f'{int(id_ann)} ({sl:.1f} m)'
            ax.annotate(id_ann, (xx_ann, yy_ann), fontsize=10, color='grey')

    ax.set_aspect(1)
    ax.set_title('Corner Reflector Locations')
    plt.gca().invert_yaxis()
    fig.savefig(outdir / f'{fileprefix}_CR_locations.png', dpi=300, bbox_inches='tight')


def calculate_ale_for_cr(point, data, outdir, fileprefix, search_window=100, oversample_factor=4):
    ybuff = search_window // 2
    xbuff = search_window // 2
    yrange = np.s_[int(point['yloc'] - ybuff) : int(point['yloc'] + ybuff)]
    xrange = np.s_[int(point['xloc'] - xbuff) : int(point['xloc'] + xbuff)]
    cropped_data = data[yrange, xrange]
    yind, xind = np.unravel_index(np.argmax(cropped_data, axis=None), cropped_data.shape)

    center_buff = 32
    yind_full = int(point['yloc'] - ybuff) + yind
    xind_full = int(point['xloc'] - xbuff) + xind
    ycenter = np.s_[int(yind_full - center_buff) : int(yind_full + center_buff)]
    xcenter = np.s_[int(xind_full - center_buff) : int(xind_full + center_buff)]
    centered_data = data[ycenter, xcenter]

    data_ovs = isce3.cal.point_target_info.oversample(centered_data, oversample_factor, baseband=True)

    yoff2 = int(data_ovs.shape[0] / 2)
    xoff2 = int(data_ovs.shape[1] / 2)
    numpix = 8
    zoom_half_size = numpix * oversample_factor
    data_ovs_zoom = data_ovs[
        yoff2 - zoom_half_size : yoff2 + zoom_half_size, xoff2 - zoom_half_size : xoff2 + zoom_half_size
    ]

    N = numpix * 2 * oversample_factor
    x = np.linspace(0, numpix * 2 * oversample_factor - 1, N)
    y = np.linspace(0, numpix * 2 * oversample_factor - 1, N)
    Xg, Yg = np.meshgrid(x, y)
    fmodel = Model(gaussfit, independent_vars=('x', 'y'))
    theta = 0.1  # deg
    x0 = numpix * oversample_factor
    y0 = numpix * oversample_factor
    sigx = 2
    sigy = 5
    A = np.max(data_ovs_zoom)
    result = fmodel.fit(data_ovs_zoom, x=Xg, y=Yg, A=A, x0=x0, y0=y0, sigma_x=sigx, sigma_y=sigy, theta=theta)
    fit = fmodel.func(Xg, Yg, **result.best_values)

    ypeak_ovs = result.best_values['y0'] + yoff2 - zoom_half_size
    ypeak_centered = ypeak_ovs / oversample_factor
    ypeak = ypeak_centered + yind_full - center_buff
    point['yloc_cr'] = ypeak

    xpeak_ovs = result.best_values['x0'] + xoff2 - zoom_half_size
    xpeak_centered = xpeak_ovs / oversample_factor
    xpeak = xpeak_centered + xind_full - center_buff
    point['xloc_cr'] = xpeak

    xreal_centered = point['xloc'] - xind_full + center_buff
    xreal_ovs = xreal_centered * oversample_factor

    yreal_centered = point['yloc'] - yind_full + center_buff
    yreal_ovs = yreal_centered * oversample_factor

    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(1, 3, figsize=(15, 7))
    ax[0].imshow(centered_data, cmap='gray', interpolation=None, origin='upper')
    ax[0].plot(xpeak_centered, ypeak_centered, 'r+', label='Return Peak')
    ax[0].plot(xreal_centered, yreal_centered, 'b+', label='CR Location')
    ax[0].legend()
    ax[0].set_title(f'Corner Reflector ID: {int(point["ID"])}')
    ax[1].imshow(data_ovs, cmap='gray', interpolation=None, origin='upper')
    ax[1].plot(xpeak_ovs, ypeak_ovs, 'r+')
    ax[1].plot(xreal_ovs, yreal_ovs, 'b+')
    ax[1].set_title(f'Oversampled Corner Reflector ID: {point["ID"]}')
    ax[2].imshow(fit, cmap='gray', interpolation=None, origin='upper')
    ax[2].plot(result.best_values['x0'], result.best_values['y0'], 'r+')
    ax[2].set_title(f'Gaussian Fit Corner Reflector ID: {int(point["ID"])}')
    [axi.axis('off') for axi in ax]
    fig.tight_layout()
    fig.savefig(outdir / f'{fileprefix}_CR_{int(point["ID"])}.png', dpi=300, bbox_inches='tight')

    return point


def cr_mean(data):
    return np.round(np.nanmean(data), 3)


def cr_spread(data):
    return np.round(np.nanstd(data) / np.sqrt(np.size(data)), 3)


def plot_ale(cr_df, azmangle, outdir, fileprefix):
    east_ale = cr_df['easting_ale']
    north_ale = cr_df['northing_ale']
    ale = cr_df['ale']
    los = np.deg2rad(90 - azmangle)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(east_ale, north_ale, s=20, c='k', alpha=0.6, marker='o')
    ax.annotate(
        'LOS',
        xytext=(np.cos(los) * 10, np.sin(los) * 10),
        xy=(0, 0),
        arrowprops=dict(edgecolor='darkblue', arrowstyle='<-'),
        color='darkblue',
    )
    ax.grid(True)
    ax.set_xlim(-15.25, 15.25)
    ax.set_ylim(-15.25, 15.25)
    ax.axhline(0, color='black')
    ax.axvline(0, color='black')
    east_metric = f'Easting: {cr_mean(east_ale)} +/- {cr_spread(east_ale)} m'
    north_metric = f'Northing: {cr_mean(north_ale)} +/- {cr_spread(north_ale)} m'
    overall_metric = f'Overall: {cr_mean(ale)} +/- {cr_spread(ale)} m'
    ax.set_title(f'{east_metric}, {north_metric}, {overall_metric}', fontsize=10)
    ax.set_xlabel('Easting Error (m)')
    ax.set_ylabel('Northing Error (m)')
    fig.suptitle('Absolute Location Error')
    plt.savefig(outdir / f'{fileprefix}_ale.png', dpi=300, bbox_inches='tight', transparent=True)


def ale(filepath, date, azmangle, outdir, fileprefix):
    outdir.mkdir(parents=True, exist_ok=True)
    ds = gdal.Open(str(filepath))
    data = ds.GetRasterBand(1).ReadAsArray()
    geotransform = ds.GetGeoTransform()
    x_start = geotransform[0] + 0.5 * geotransform[1]
    y_start = geotransform[3] + 0.5 * geotransform[5]
    x_end = x_start + geotransform[1] * ds.RasterXSize
    y_end = y_start + geotransform[5] * ds.RasterYSize
    bounds = (x_start, y_start, x_end, y_end)
    bounds = box(*bounds)
    x_spacing = geotransform[1]
    y_spacing = geotransform[5]

    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjectionRef())
    epsg = int(srs.GetAuthorityCode(None))
    cr_df = get_cr_df(bounds, epsg, date, outdir)
    cr_df = add_image_location(cr_df, epsg, x_start, y_start, x_spacing, y_spacing, bounds)
    cr_df = filter_valid_data(cr_df, data)
    cr_df = filter_orientation(cr_df, azmangle)
    cr_df = cr_df.assign(yloc_cr=np.nan, xloc_cr=np.nan)
    plot_crs_on_image(cr_df, data, outdir, fileprefix)
    for idx, cr in cr_df.iterrows():
        cr = calculate_ale_for_cr(cr, data, outdir, fileprefix)
        cr_df.iloc[idx] = cr

    cr_df['easting_ale'] = (cr_df['xloc_cr'] - cr_df['xloc_floats']) * x_spacing
    cr_df['northing_ale'] = (cr_df['yloc_cr'] - cr_df['yloc_floats']) * y_spacing
    cr_df['ale'] = np.sqrt(cr_df['northing_ale'] ** 2 + cr_df['easting_ale'] ** 2)
    cr_df.to_csv(outdir / (fileprefix + '_ale.csv'), index=False)
    plot_ale(cr_df, azmangle, outdir, fileprefix)


def main():
    """Example:
    python ale.py rtc.tif 2024-12-03 DESC out --fileprefix rtc
    """
    parser = ArgumentParser(description='Absolute Location Error Estimation.')
    parser.add_argument('filepath', type=str, help='Path to the file to be processed')
    parser.add_argument('date', type=str, help='Date of the image collection (YYYY-MM-DD)')
    parser.add_argument('azmangle', type=int, help='Azimuth angle of the image (clockwise from North in degrees)')
    parser.add_argument('outdir', type=str, help='Directory to save the results')
    parser.add_argument(
        '--fileprefix', type=str, help='Prefix for the output filenames (default: input filename)', default=None
    )
    args = parser.parse_args()
    args.filepath = Path(args.filepath)
    args.date = datetime.strptime(args.date, '%Y-%m-%d')
    assert 0 <= args.azmangle <= 360, f'Azimuth angle {args.azmangle} is out of range [0, 360].'
    args.outdir = Path(args.outdir)
    if args.fileprefix is None:
        args.fileprefix = Path(args.filepath).stem
    assert args.filepath.exists(), f'File {args.filepath} does not exist.'

    ale(args.filepath, args.date, args.azmangle, args.outdir, args.fileprefix)


if __name__ == '__main__':
    main()
