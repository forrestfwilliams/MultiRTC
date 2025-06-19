import numpy as np
import pandas as pd
import pyproj
import requests
from pyproj import Transformer, transform
from shapely.geometry import Point, box


def filter_orientation(cr_df, azm_angle):
    looking_east = azm_angle < 180
    if looking_east:
        cr_df = cr_df[(cr_df['azm'] < 200) & (cr_df['azm'] > 20)].reset_index(drop=True)
    else:
        cr_df = cr_df[cr_df['azm'] > 340].reset_index(drop=True)
    return cr_df


def get_cr_df(bounds, epsg, date, azm_angle, outdir):
    rosamond_bounds = box(*[-124.409591, 32.534156, -114.131211, 42.009518])
    if epsg != 4326:
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
    not_in_bounds = []
    for idx, row in cr_df.iterrows():
        point = Point(row['lon'], row['lat'])
        if not bounds.contains(point):
            not_in_bounds.append(idx)
    cr_df = cr_df.drop(not_in_bounds).reset_index(drop=True)
    cr_df = filter_orientation(cr_df, azm_angle)
    return cr_df


def add_geo_image_location(cr_df, epsg, x_start, y_start, x_spacing, y_spacing, bounds):
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


def add_rdr_image_location(slc, cr_df, search_radius):
    blank = [np.nan] * cr_df.shape[0]
    cr_df = cr_df.assign(xloc=blank, yloc=blank)
    llh2ecef = pyproj.Transformer.from_crs('EPSG:4979', 'EPSG:4978', always_xy=True)
    no_peak = []
    for idx, row in cr_df.iterrows():
        x, y, z = llh2ecef.transform(row['lon'], row['lat'], slc.scp_hae)
        row_guess, col_guess = slc.geo2rowcol(np.array([[x, y, z]]))[0]
        row_guess, col_guess = int(round(row_guess)), int(round(col_guess))
        row_range = (row_guess - search_radius, row_guess + search_radius)
        col_range = (col_guess - search_radius, col_guess + search_radius)
        data = slc.load_data(row_range, col_range)
        in_db = 10 * np.log10(np.abs(data))
        row_peak, col_peak = np.unravel_index(np.argmax(in_db, axis=None), in_db.shape)
        if in_db[row_peak, col_peak] < 20:
            no_peak.append(idx)
            print(f'No peak found for CR {int(row["ID"])}')
            continue
        row['yloc'] = row_range[0] + row_peak
        row['xloc'] = col_range[0] + col_peak
        cr_df.iloc[idx] = row

    cr_df = cr_df.drop(no_peak).reset_index(drop=True)
    return cr_df
