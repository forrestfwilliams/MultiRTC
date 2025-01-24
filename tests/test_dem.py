from shapely.geometry import box

from multirtc import dem


def test_get_granule_url():
    test_url = 'https://nisar.asf.earthdatacloud.nasa.gov/STATIC/DEM/v1.1/EPSG4326/S10_W020/DEM_S01_00_W001_00.tif'
    url = dem.get_dem_granule_url(-1, -1)
    assert url == test_url


def test_get_latlon_pairs():
    polygon = box(-1, -1, 1, 1)
    latlon_pairs = dem.get_latlon_pairs(polygon)
    assert latlon_pairs == [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
