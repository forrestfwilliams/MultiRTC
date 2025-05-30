import isce3
import numpy as np
import pyproj
from shapely.geometry import Polygon


ECEF = pyproj.CRS(4978)
LLA = pyproj.CRS(4979)
ECEF2LLA = pyproj.Transformer.from_crs(ECEF, LLA, always_xy=True)


def get_point_epsg(lat, lon):
    if (lon >= 180.0) or (lon <= -180.0):
        lon = (lon + 180.0) % 360.0 - 180.0
    if lat >= 75.0:
        epsg = 3413
    elif lat <= -75.0:
        epsg = 3031
    elif lat > 0:
        epsg = 32601 + int(np.round((lon + 177) / 6.0))
    elif lat < 0:
        epsg = 32701 + int(np.round((lon + 177) / 6.0))
    else:
        raise ValueError(f'Could not determine EPSG for {lon}, {lat}')
    assert 1024 <= epsg <= 32767, 'Computed EPSG is out of range'
    return epsg


def snap_coord(val, snap, round_func):
    """
    Returns the snapped values of the input value

    Parameters
    -----------
    val : float
        Input value to snap
    snap : float
        Snapping step
    round_func : function pointer
        A function used to round `val` i.e. round, ceil, floor

    Return:
    --------
    snapped_value : float
        snapped value of `var` by `snap`

    """
    snapped_value = round_func(float(val) / snap) * snap
    return snapped_value


def grid_size(stop, start, sz):
    """
    get grid dim based on start, end, and grid size inputs
    """
    assert None not in [stop, start, sz], 'Invalid input values'
    return int(np.round(np.abs((stop - start) / sz)))


def snap_geogrid(geogrid, x_snap, y_snap):
    """
    Snap geogrid based on user-defined snapping values

    Parameters
    ----------
    geogrid: isce3.product.GeoGridParameters
        ISCE3 object definining the geogrid
    x_snap: float
        Snap value along X-direction
    y_snap: float
        Snap value along Y-direction

    Returns
    -------
    geogrid: isce3.product.GeoGridParameters
        ISCE3 object containing the snapped geogrid
    """
    xmax = geogrid.start_x + geogrid.width * geogrid.spacing_x
    ymin = geogrid.start_y + geogrid.length * geogrid.spacing_y

    geogrid.start_x = snap_coord(geogrid.start_x, x_snap, np.floor)
    end_x = snap_coord(xmax, x_snap, np.ceil)
    geogrid.width = grid_size(end_x, geogrid.start_x, geogrid.spacing_x)

    geogrid.start_y = snap_coord(geogrid.start_y, y_snap, np.ceil)
    end_y = snap_coord(ymin, y_snap, np.floor)
    geogrid.length = grid_size(end_y, geogrid.start_y, geogrid.spacing_y)
    return geogrid


def get_geogrid_poly(geogrid):
    new_maxx = geogrid.start_x + (geogrid.width * geogrid.spacing_x)
    new_miny = geogrid.start_y + (geogrid.length * geogrid.spacing_y)
    points = [
        [geogrid.start_x, geogrid.start_y],
        [geogrid.start_x, new_miny],
        [new_maxx, new_miny],
        [new_maxx, geogrid.start_y],
    ]
    poly = Polygon(points)
    return poly


def generate_geogrids(slc_obj, resolution: int, epsg: int = None, rda: bool = True):
    """
    Compute the slc geogrid
    """
    x_spacing = resolution
    y_spacing = -1 * np.abs(resolution)
    if epsg is None:
        epsg = get_point_epsg(slc_obj.center.y, slc_obj.center.x)

    # if rda:
    #     radar_grid = slc_obj.as_isce3_radargrid()
    # else:
    #     geogrid = slc_obj.get_geogrid(x_spacing, y_spacing)
    geogrid = isce3.product.bbox_to_geogrid(
        slc_obj.radar_grid, slc_obj.orbit, slc_obj.doppler_centroid_grid, x_spacing, y_spacing, epsg
    )
    geogrid_snapped = snap_geogrid(geogrid, geogrid.spacing_x, geogrid.spacing_y)
    return geogrid_snapped
