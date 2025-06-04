import isce3
import numpy as np
from shapely.geometry import Polygon


def get_point_epsg(lat: float, lon: float) -> int:
    """Determine the best EPSG code for a given latitude and longitude.
    Returns the local UTM zone for latitudes between +/-75 degrees and
    polar/antartic stereographic for latitudes outside that range.

    Args:
        lat: Latitude in degrees.
        lon: Longitude in degrees.

    Returns:
        EPSG code for the specified latitude and longitude.
    """
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
    assert (32600 <= epsg <= 32761) or epsg in [3031, 3413], 'Computed EPSG is out of range'
    return epsg


def snap_coord(val: float, snap: float, round_func: callable) -> float:
    """
    Returns the snapped version of the input value

    Args:
        val : value to snap
        snap : snapping step
        round_func : function pointer to round, ceil, or floor

    Returns:
        snapped value of `val` by `snap`
    """
    snapped_value = round_func(float(val) / snap) * snap
    return snapped_value


def grid_size(stop: float, start: float, size: float):
    """
    Get number of grid points based on start, end, and grid size inputs

    Args:
        stop: End value of grid
        start: Start value of grid
        size: Grid size in same units as start and stop

    Returns:
        Number of grid points between start and stop
    """
    return int(np.round(np.abs((stop - start) / size)))


def snap_geogrid(
    geogrid: isce3.product.GeoGridParameters, x_snap: float, y_snap: float
) -> isce3.product.GeoGridParameters:
    """
    Snap geogrid based on user-defined snapping values

    Args:
        geogrid: ISCE3 object definining the geogrid
        x_snap: Snap value along X-direction
        y_snap: Snap value along Y-direction

    Returns:
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


def get_geogrid_poly(geogrid: isce3.product.GeoGridParameters) -> Polygon:
    """
    Create a polygon from a geogrid object

    Args:
        geogrid: ISCE3 object defining the geogrid

    Returns:
        Shapely Polygon representing the geogrid area
    """
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


def generate_geogrids(slc, spacing_meters: int, epsg: int) -> isce3.product.GeoGridParameters:
    """Compute a geogrid based on the radar grid of the SLC and the specified spacing.

    Args:
        slc: Slc-derived object containing radar grid, orbit, and doppler centroid grid.
        spacing_meters: Spacing in meters for the geogrid.
        epsg: EPSG code for the coordinate reference system.

    Returns:
        A geogrid object with the specified spacing.
    """
    x_spacing = spacing_meters
    y_spacing = -1 * np.abs(spacing_meters)
    geogrid = isce3.product.bbox_to_geogrid(
        slc.radar_grid, slc.orbit, slc.doppler_centroid_grid, x_spacing, y_spacing, epsg
    )
    geogrid_snapped = snap_geogrid(geogrid, geogrid.spacing_x, geogrid.spacing_y)
    return geogrid_snapped
