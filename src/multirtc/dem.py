"""Prepare an external DEM for use in MultiRTC"""

from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path
from shutil import copyfile
from tempfile import NamedTemporaryFile

import numpy as np
import shapely
from osgeo import gdal, osr
from shapely.geometry import LinearRing, Polygon, box

from multirtc.fetch import download_file


gdal.UseExceptions()
URL = 'https://nisar.asf.earthdatacloud.nasa.gov/STATIC/DEM/v1.1/EPSG4326'
EGM2008_GEOID = {
    'WORLD': [
        '/vsicurl/https://asf-dem-west.s3.amazonaws.com/GEOID/us_nga_egm2008_1.tif',
        box(-180.0083333, -90.0083333, 180.0083333, 90.0083333),
    ]
}
NAD88 = {
    'CONUS': [
        '/vsicurl/https://asf-dem-west.s3.amazonaws.com/GEOID/us_noaa_g2012bu0_wgs84.tif',
        box(-130.0083313, 23.9916667, -59.9920366, 58.0083351),
    ],
    'AK_EAST': [
        '/vsicurl/https://asf-dem-west.s3.amazonaws.com/GEOID/us_noaa_g2012ba0_wgs84_east.tif',
        box(171.9916687, 48.9916583, 234.008, 72.008),
    ],
    'AK_WEST': [
        '/vsicurl/https://asf-dem-west.s3.amazonaws.com/GEOID/us_noaa_g2012ba0_wgs84_west.tif',
        box(-188.008, 48.992, -125.9915921, 72.0083313),
    ],
}
VALID_DATUMS = ['WGS84', 'EGM2008', 'NAD88']


def get_bbox_from_info(info: dict) -> Polygon:
    minx = info['cornerCoordinates']['lowerLeft'][0]
    miny = info['cornerCoordinates']['lowerLeft'][1]
    maxx = info['cornerCoordinates']['upperRight'][0]
    maxy = info['cornerCoordinates']['upperRight'][1]
    return box(minx, miny, maxx, maxy)


def validate_dem(dem_path: Path, footprint: Polygon) -> None:
    """Validate that the DEM file is in EPSG:4326 and contains the given footprint.

    Args:
        dem_path: Path to the DEM file.
        footprint: Polygon representing the area of interest.

    Raises:
        ValueError: If the DEM does not cover the footprint.
    """
    info = gdal.Info(str(dem_path), format='json')

    srs = osr.SpatialReference()
    srs.ImportFromWkt(info['coordinateSystem']['wkt'])
    assert int(srs.GetAttrValue('AUTHORITY', 1)) == 4326, f'DEM file {dem_path} is not in EPSG:4326 projection.'

    dem_extent = get_bbox_from_info(info)
    if not dem_extent.contains(footprint):
        dem_bound_str = ', '.join([str(round(x, 3)) for x in dem_extent.bounds])
        footprint_bound_str = ', '.join([str(round(x, 3)) for x in footprint.bounds])
        raise ValueError(
            f'DEM does not fully cover the footprint: ({dem_bound_str}) for DEM, vs ({footprint_bound_str})'
        )


def check_antimeridean(poly: Polygon) -> list[Polygon]:
    """Check if the polygon crosses the antimeridian and split the polygon if it does.

    Args:
        poly: Polygon object to check for antimeridian crossing.

    Returns:
        List of Polygon objects, split if necessary.
    """
    x_min, _, x_max, _ = poly.bounds

    # Check anitmeridean crossing
    if (x_max - x_min > 180.0) or (x_min <= 180.0 <= x_max):
        dateline = shapely.wkt.loads('LINESTRING( 180.0 -90.0, 180.0 90.0)')

        # build new polygon with all longitudes between 0 and 360
        x, y = poly.exterior.coords.xy
        new_x = (k + (k <= 0.0) * 360 for k in x)
        new_ring = LinearRing(zip(new_x, y))

        # Split input polygon
        # (https://gis.stackexchange.com/questions/232771/splitting-polygon-by-linestring-in-geodjango_)
        merged_lines = shapely.ops.linemerge([dateline, new_ring])
        border_lines = shapely.ops.unary_union(merged_lines)
        decomp = shapely.ops.polygonize(border_lines)

        polys = list(decomp)

        for polygon_count in range(len(polys)):
            x, y = polys[polygon_count].exterior.coords.xy
            # if there are no longitude values above 180, continue
            if not any([k > 180 for k in x]):
                continue

            # otherwise, wrap longitude values down by 360 degrees
            x_wrapped_minus_360 = np.asarray(x) - 360
            polys[polygon_count] = Polygon(zip(x_wrapped_minus_360, y))

    else:
        # If dateline is not crossed, treat input poly as list
        polys = [poly]

    return polys


def get_dem_granule_url(lat: int, lon: int) -> str:
    """Generate the URL for the OPERA DEM granule based on latitude and longitude.

    Args:
        lat: Latitude in degrees.
        lon: Longitude in degrees.

    Returns:
        URL string for the DEM granule.
    """
    lat_tens = np.floor_divide(lat, 10) * 10
    lat_cardinal = 'S' if lat_tens < 0 else 'N'

    lon_tens = np.floor_divide(lon, 20) * 20
    lon_cardinal = 'W' if lon_tens < 0 else 'E'

    prefix = f'{lat_cardinal}{np.abs(lat_tens):02d}_{lon_cardinal}{np.abs(lon_tens):03d}'
    filename = f'DEM_{lat_cardinal}{np.abs(lat):02d}_00_{lon_cardinal}{np.abs(lon):03d}_00.tif'
    file_url = f'{URL}/{prefix}/{filename}'
    return file_url


def get_latlon_pairs(polygon: Polygon) -> list[tuple[float, float]]:
    """Get latitude and longitude pairs for the bounding box of a polygon.

    Args:
        polygon: Polygon object representing the area of interest.

    Returns:
        List of tuples containing latitude and longitude pairs for each point of the bounding box.
    """
    minx, miny, maxx, maxy = polygon.bounds
    lats = np.arange(np.floor(miny), np.floor(maxy) + 1).astype(int)
    lons = np.arange(np.floor(minx), np.floor(maxx) + 1).astype(int)
    return list(product(lats, lons))


def download_opera_dem_for_footprint(output_path: Path, footprint: Polygon, buffer: float = 0.2) -> None:
    """
    Download the OPERA DEM for a given footprint and save it to the specified output path.

    Args:
        output_path: Path where the DEM will be saved.
        footprint: Polygon representing the area of interest.
        buffer: Buffer distance in degrees to extend the footprint.
    """
    output_dir = output_path.parent
    if output_path.exists():
        return output_path

    footprint = box(*footprint.buffer(buffer).bounds)
    footprints = check_antimeridean(footprint)
    latlon_pairs = []
    for footprint in footprints:
        latlon_pairs += get_latlon_pairs(footprint)
    urls = [get_dem_granule_url(lat, lon) for lat, lon in latlon_pairs]

    with ThreadPoolExecutor(max_workers=4) as executor:
        executor.map(lambda url: download_file(url, str(output_dir)), urls)

    # [download_file(url, str(output_dir)) for url in urls]
    vrt_filepath = output_dir / 'dem.vrt'
    input_files = [str(output_dir / Path(url).name) for url in urls]
    gdal.BuildVRT(str(output_dir / 'dem.vrt'), input_files)
    ds = gdal.Open(str(vrt_filepath), gdal.GA_ReadOnly)
    gdal.Translate(str(output_path), ds, format='GTiff')

    ds = None
    [Path(f).unlink() for f in input_files + [vrt_filepath]]


def get_correction_geoid(bbox: Polygon, input_datum: str) -> str:
    """Get the path to the geoid correction file based on the bounding box and input datum."""
    if input_datum.upper() == 'EGM2008':
        return EGM2008_GEOID['WORLD'][0]

    if input_datum.upper() == 'NAD88':
        for region, (correction_path, region_bbox) in NAD88.items():
            if bbox.intersects(region_bbox):
                return correction_path

    raise ValueError(f'No suitable geoid correction found for datum {input_datum} and bounding box {bbox}.')


def convert_to_height_above_ellipsoid(dem_file: Path, input_datum) -> None:
    assert input_datum in VALID_DATUMS, f'Input datum must be one of {VALID_DATUMS}, got {input_datum}.'
    if input_datum.upper() == 'WGS84':
        return
    dem_info = gdal.Info(str(dem_file), format='json')
    dem_bbox = get_bbox_from_info(dem_info)
    correction_path = get_correction_geoid(dem_bbox, input_datum)
    with NamedTemporaryFile() as geoid_file:
        gdal.Warp(
            geoid_file.name,
            correction_path,
            dstSRS=dem_info['coordinateSystem']['wkt'],
            outputBounds=dem_bbox.bounds,
            width=dem_info['size'][0],
            height=dem_info['size'][1],
            resampleAlg='cubic',
            multithread=True,
            format='GTiff',
        )
        geoid_ds = gdal.Open(geoid_file.name)
        geoid_data = geoid_ds.GetRasterBand(1).ReadAsArray()
        del geoid_ds

        dem_ds = gdal.Open(str(dem_file), gdal.GA_Update)
        dem_data = dem_ds.GetRasterBand(1).ReadAsArray()
        nan_value = dem_ds.GetRasterBand(1).GetNoDataValue()
        if nan_value is not None:
            nan_mask = dem_data == nan_value
            dem_data += geoid_data
            dem_data[nan_mask] = nan_value
        else:
            dem_data += geoid_data
        dem_ds.GetRasterBand(1).WriteArray(dem_data)
        dem_ds.FlushCache()
        del dem_ds


def set_nodata_value(dem_file: Path, nodata_value: float) -> None:
    """Set the NoData value for the DEM file.

    Args:
        dem_file: Path to the DEM file.
        nodata_value: Value to set as NoData.
    """
    dem_ds = gdal.Open(str(dem_file), gdal.GA_Update)
    band = dem_ds.GetRasterBand(1)
    if band.GetNoDataValue() is not None:
        data = band.ReadAsArray()
        data[data == band.GetNoDataValue()] = nodata_value
        band.WriteArray(data)
    band.SetNoDataValue(nodata_value)
    band.FlushCache()
    del dem_ds


def prep_dem(input_path: Path, output_path: Path, input_datum: str) -> None:
    """Prepare the DEM for processing by reprojecting it to EPSG:4326 and (optionally) converting it to height above ellipsoid.

    Args:
        dem_path: Path to the DEM file.
        input_datum: Datum of the input DEM, either 'WGS84', 'EGM2008', 'NAD88'.
        output_path: Path where the prepared DEM will be saved.
    """
    assert input_datum in VALID_DATUMS, f'Input datum must be one of {VALID_DATUMS}, got {input_datum}.'
    copyfile(input_path, output_path)
    info = gdal.Info(str(input_path), format='json')
    srs = osr.SpatialReference()
    srs.ImportFromWkt(info['coordinateSystem']['wkt'])
    if int(srs.GetAttrValue('AUTHORITY', 1)) != 4326:
        gdal.Warp(
            str(output_path),
            str(output_path),
            dstSRS='EPSG:4326',
            resampleAlg='cubic',
            multithread=True,
        )
    convert_to_height_above_ellipsoid(output_path, input_datum.upper())
    set_nodata_value(output_path, 0)


def create_parser(parser):
    parser.add_argument('input_dem', help='Path to the input DEM file to prepare for processing.')
    parser.add_argument('output_dem', help='Path where the prepared DEM will be saved.')
    parser.add_argument(
        'input_datum', choices=VALID_DATUMS, help='Datum of the input DEM (supported: WGS84, EGM2008, NAD88).'
    )
    return parser


def run(args):
    args.input_dem = Path(args.input_dem)
    assert args.input_dem.exists(), f'DEM file {args.input_dem} does not exist.'
    args.output_dem = Path(args.output_dem)
    prep_dem(args.input_dem, args.output_dem, args.input_datum)
