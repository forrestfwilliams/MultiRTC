from pathlib import Path

import numpy as np
from osgeo import gdal
from pyproj import Transformer


gdal.UseExceptions()


def write_epsg4326_geotiff(output_path: Path, data: np.ndarray, transform: tuple) -> None:
    """Write a numpy array to a GeoTIFF file.

    Args:
        output_path: Path to the output GeoTIFF file.
        data: 2D numpy array containing the data to write.
        transform: GeoTransform tuple for the raster.
    """
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(str(output_path), data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform(transform)
    ds.SetProjection('EPSG:4326')  # WGS84
    ds.GetRasterBand(1).WriteArray(data)
    ds.FlushCache()
    ds = None


def convert_from_nad83_to_wgs84(navd83: str) -> None:
    """Convert NAD83(2011) ellipsoidal heights to WGS84 ellipsoidal heights.

    Args:
        navd83: Path to the NAD83(2011) geoid model.
    """
    ds = gdal.Open(navd83, gdal.GA_ReadOnly)
    nad83_heights = ds.GetRasterBand(1).ReadAsArray()
    transform = ds.GetGeoTransform()
    lons = np.zeros_like(nad83_heights)
    lats = np.zeros_like(nad83_heights)
    wgs84_heights = np.zeros_like(nad83_heights)
    ds = None

    # NAD83(2011) ellipsoidal heights → WGS84 ellipsoidal heights
    transformer = Transformer.from_crs('EPSG:6319', 'EPSG:4979', always_xy=True)
    for i in range(nad83_heights.shape[0]):
        for j in range(nad83_heights.shape[1]):
            lon, lat = transform[0] + j * transform[1], transform[3] + i * transform[5]
            lon_wgs84, lat_wgs84, h_wgs84 = transformer.transform(lon, lat, nad83_heights[i, j])
            lons[i, j] = lon_wgs84
            lats[i, j] = lat_wgs84
            wgs84_heights[i, j] = h_wgs84

    wgs84_lon_start = lons[0, 0]
    wgs84_lon_step = np.diff(lons, axis=1).mean()
    wgs84_lat_start = lats[0, 0]
    wgs84_lat_step = np.diff(lats, axis=0).mean()
    trans = [wgs84_lon_start, wgs84_lon_step, 0, wgs84_lat_start, 0, wgs84_lat_step]
    if trans[0] > 180:
        trans[0] -= 360
    if (trans[0] + (wgs84_lon_step * wgs84_heights.shape[1])) > 180:
        write_epsg4326_geotiff(f'{Path(navd83).stem}_wgs84_east.tif', wgs84_heights, trans)
        west_trans = [trans[0] - 360, trans[1], trans[2], trans[3], trans[4], trans[5]]
        write_epsg4326_geotiff(f'{Path(navd83).stem}_wgs84_west.tif', wgs84_heights, west_trans)
    else:
        write_epsg4326_geotiff(f'{Path(navd83).stem}_wgs84.tif', wgs84_heights, trans)


convert_from_nad83_to_wgs84(
    '/vsicurl/https://github.com/OSGeo/PROJ-data/raw/refs/heads/master/us_noaa/us_noaa_g2012ba0.tif'
)
convert_from_nad83_to_wgs84(
    '/vsicurl/https://github.com/OSGeo/PROJ-data/raw/refs/heads/master/us_noaa/us_noaa_g2012bu0.tif'
)
