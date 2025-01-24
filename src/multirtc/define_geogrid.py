import isce3
import numpy as np


def get_point_epsg(lat, lon):
    """
    Get EPSG code based on latitude and longitude
    coordinates of a point

    Parameters
    ----------
    lat: float
        Latitude coordinate of the point
    lon: float
        Longitude coordinate of the point

    Returns
    -------
    epsg: int
        UTM zone, Polar stereographic (North / South)
    """
    # "Warp" longitude value into [-180.0, 180.0]
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


def generate_geogrids(burst, x_spacing, y_spacing, x_snap, y_snap, output_epsg=None):
    """
    Compute frame and bursts geogrids

    Parameters
    ----------
    bursts: list[s1reader.s1_burst_slc.Sentinel1BurstSlc]
        List of S-1 burst SLC objects
    opts: RtcOptions

    Returns
    -------
    geogrid_all_snapped: isce3.product.GeoGridParameters
        Mosaic geogrid
    geogrids_dict: dict
        Dict containing bursts' geogrids indexed by burst_id
    """
    if output_epsg is None:
        output_epsg = get_point_epsg(burst.center.y, burst.center.x)
    epsg_bursts = output_epsg
    y_spacing_negative = -1 * np.abs(y_spacing)

    radar_grid = burst.as_isce3_radargrid()
    geogrid = isce3.product.bbox_to_geogrid(
        radar_grid, burst.orbit, isce3.core.LUT2d(), x_spacing, y_spacing_negative, epsg_bursts
    )
    geogrid_snapped = snap_geogrid(geogrid, x_snap, y_snap)
    return geogrid_snapped
