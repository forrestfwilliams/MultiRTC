import logging

import isce3
import numpy as np
from osgeo import gdal


logger = logging.getLogger('rtc_s1')


def compute_correction_lut(
    burst,
    dem_raster,
    scratch_path,
    rg_step_meters,
    az_step_meters,
    apply_bistatic_delay_correction,
    apply_static_tropospheric_delay_correction,
):
    """
    Compute lookup table for geolocation correction.
    Applied corrections are: bistatic delay (azimuth),
                             static troposphere delay (range)

    Parameters
    ----------
    burst_in: Sentinel1BurstSlc
        Input burst SLC
    dem_raster: isce3.io.raster
        DEM to run rdr2geo
    scratch_path: str
        Scratch path where the radargrid rasters will be saved
    rg_step_meters: float
        LUT spacing in slant range. Unit: meters
    az_step_meters: float
        LUT spacing in azimth direction. Unit: meters
    apply_bistatic_delay_correction: bool
        Flag to indicate whether the bistatic delay correciton should be applied
    apply_static_tropospheric_delay_correction: bool
        Flag to indicate whether the static tropospheric delay correction should be
        applied

    Returns
    -------
    rg_lut, az_lut: isce3.core.LUT2d
        LUT2d for geolocation correction in slant range and azimuth direction
    """

    rg_lut = None
    az_lut = None

    # approximate conversion of az_step_meters from meters to seconds
    numrow_orbit = burst.orbit.position.shape[0]
    vel_mid = burst.orbit.velocity[numrow_orbit // 2, :]
    spd_mid = np.linalg.norm(vel_mid)
    pos_mid = burst.orbit.position[numrow_orbit // 2, :]
    alt_mid = np.linalg.norm(pos_mid)

    r = 6371000.0  # geometric mean of WGS84 ellipsoid

    az_step_sec = (az_step_meters * alt_mid) / (spd_mid * r)
    # Bistatic - azimuth direction
    bistatic_delay = burst.bistatic_delay(range_step=rg_step_meters, az_step=az_step_sec)

    if apply_bistatic_delay_correction:
        az_lut = isce3.core.LUT2d(
            bistatic_delay.x_start,
            bistatic_delay.y_start,
            bistatic_delay.x_spacing,
            bistatic_delay.y_spacing,
            -bistatic_delay.data,
        )

    if not apply_static_tropospheric_delay_correction:
        return rg_lut, az_lut

    # Calculate rdr2geo rasters
    epsg = dem_raster.get_epsg()
    proj = isce3.core.make_projection(epsg)
    ellipsoid = proj.ellipsoid

    rdr_grid = burst.as_isce3_radargrid(az_step=az_step_sec, rg_step=rg_step_meters)

    grid_doppler = isce3.core.LUT2d()

    # Initialize the rdr2geo object
    rdr2geo_obj = isce3.geometry.Rdr2Geo(rdr_grid, burst.orbit, ellipsoid, grid_doppler, threshold=1.0e-8)

    # Get the rdr2geo raster needed for SET computation
    topo_output = {
        f'{scratch_path}/height.rdr': gdal.GDT_Float32,
        f'{scratch_path}/incidence_angle.rdr': gdal.GDT_Float32,
    }

    raster_list = []
    for fname, dtype in topo_output.items():
        topo_output_raster = isce3.io.Raster(fname, rdr_grid.width, rdr_grid.length, 1, dtype, 'ENVI')
        raster_list.append(topo_output_raster)

    height_raster, incidence_raster = raster_list

    rdr2geo_obj.topo(
        dem_raster, x_raster=None, y_raster=None, height_raster=height_raster, incidence_angle_raster=incidence_raster
    )

    height_raster.close_dataset()
    incidence_raster.close_dataset()

    # Load height and incidence angle layers
    height_arr = gdal.Open(f'{scratch_path}/height.rdr', gdal.GA_ReadOnly).ReadAsArray()
    incidence_angle_arr = gdal.Open(f'{scratch_path}/incidence_angle.rdr', gdal.GA_ReadOnly).ReadAsArray()

    # static troposphere delay - range direction
    # reference:
    # Breit et al., 2010, TerraSAR-X SAR Processing and Products,
    # IEEE Transactions on Geoscience and Remote Sensing, 48(2), 727-740.
    # DOI: 10.1109/TGRS.2009.2035497
    zenith_path_delay = 2.3
    reference_height = 6000.0
    tropo = zenith_path_delay / np.cos(np.deg2rad(incidence_angle_arr)) * np.exp(-1 * height_arr / reference_height)

    # Prepare the computation results into LUT2d
    rg_lut = isce3.core.LUT2d(
        bistatic_delay.x_start, bistatic_delay.y_start, bistatic_delay.x_spacing, bistatic_delay.y_spacing, tropo
    )

    return rg_lut, az_lut
