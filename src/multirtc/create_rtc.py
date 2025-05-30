import itertools
import logging
import os
import time

import isce3
import numpy as np
import pyproj
from osgeo import gdal
from scipy import ndimage
from tqdm import tqdm

from multirtc.prep_burst import S1BurstSlc
from multirtc.s1burst_corrections import compute_correction_lut
from multirtc.sicd import SicdSlc


logger = logging.getLogger('rtc_s1')

LAYER_NAME_LAYOVER_SHADOW_MASK = 'mask'
LAYER_NAME_RTC_ANF_GAMMA0_TO_SIGMA0 = 'rtc_anf_gamma0_to_sigma0'
LAYER_NAME_NUMBER_OF_LOOKS = 'number_of_looks'
LAYER_NAME_INCIDENCE_ANGLE = 'incidence_angle'
LAYER_NAME_LOCAL_INCIDENCE_ANGLE = 'local_incidence_angle'
LAYER_NAME_PROJECTION_ANGLE = 'projection_angle'
LAYER_NAME_RTC_ANF_PROJECTION_ANGLE = 'rtc_anf_projection_angle'
LAYER_NAME_RANGE_SLOPE = 'range_slope'
LAYER_NAME_DEM = 'interpolated_dem'


def compute_layover_shadow_mask(
    radar_grid: isce3.product.RadarGridParameters,
    orbit: isce3.core.Orbit,
    geogrid_in: isce3.product.GeoGridParameters,
    dem_raster: isce3.io.Raster,
    filename_out: str,
    output_raster_format: str,
    scratch_dir: str,
    shadow_dilation_size: int,
    threshold_rdr2geo: float = 1.0e-7,
    numiter_rdr2geo: int = 25,
    extraiter_rdr2geo: int = 10,
    lines_per_block_rdr2geo: int = 1000,
    threshold_geo2rdr: float = 1.0e-7,
    numiter_geo2rdr: int = 25,
    memory_mode: isce3.core.GeocodeMemoryMode = None,
    geocode_options=None,
    doppler=None,
):
    """
    Compute the layover/shadow mask and geocode it

    Parameters
    -----------
    radar_grid: isce3.product.RadarGridParameters
        Radar grid
    orbit: isce3.core.Orbit
        Orbit defining radar motion on input path
    geogrid_in: isce3.product.GeoGridParameters
        Geogrid to geocode the layover/shadow mask in radar grid
    geogrid_in: isce3.product.GeoGridParameters
        Geogrid to geocode the layover/shadow mask in radar grid
    dem_raster: isce3.io.Raster
        DEM raster
    filename_out: str
        Path to the geocoded layover/shadow mask
    output_raster_format: str
        File format of the layover/shadow mask
    scratch_dir: str
        Temporary Directory
    shadow_dilation_size: int
        Layover/shadow mask dilation size of shadow pixels
    threshold_rdr2geo: float
        Iteration threshold for rdr2geo
    numiter_rdr2geo: int
        Number of max. iteration for rdr2geo object
    extraiter_rdr2geo: int
        Extra number of iteration for rdr2geo object
    lines_per_block_rdr2geo: int
        Lines per block for rdr2geo
    threshold_geo2rdr: float
        Iteration threshold for geo2rdr
    numiter_geo2rdr: int
        Number of max. iteration for geo2rdr object
    memory_mode: isce3.core.GeocodeMemoryMode
        Geocoding memory mode
    geocode_options: dict
        Keyword arguments to be passed to the geocode() function
        when map projection the layover/shadow mask

    Returns
    -------
    slantrange_layover_shadow_mask_raster: isce3.io.Raster
        Layover/shadow-mask ISCE3 raster object in radar coordinates
    """
    if doppler is None:
        doppler = isce3.core.LUT2d()

    # Run topo to get layover/shadow mask
    ellipsoid = isce3.core.Ellipsoid()
    grid_doppler = doppler
    rdr2geo_obj = isce3.geometry.Rdr2Geo(
        radar_grid,
        orbit,
        ellipsoid,
        grid_doppler,
        threshold=threshold_rdr2geo,
        numiter=numiter_rdr2geo,
        extraiter=extraiter_rdr2geo,
        lines_per_block=lines_per_block_rdr2geo,
    )

    if shadow_dilation_size > 0:
        path_layover_shadow_mask_file = os.path.join(scratch_dir, 'layover_shadow_mask_slant_range.tif')
        slantrange_layover_shadow_mask_raster = isce3.io.Raster(
            path_layover_shadow_mask_file, radar_grid.width, radar_grid.length, 1, gdal.GDT_Byte, 'GTiff'
        )
    else:
        slantrange_layover_shadow_mask_raster = isce3.io.Raster(
            'layover_shadow_mask', radar_grid.width, radar_grid.length, 1, gdal.GDT_Byte, 'MEM'
        )

    rdr2geo_obj.topo(dem_raster, layover_shadow_raster=slantrange_layover_shadow_mask_raster)

    if shadow_dilation_size > 1:
        """
        constants from ISCE3:
            SHADOW_VALUE = 1;
            LAYOVER_VALUE = 2;
            LAYOVER_AND_SHADOW_VALUE = 3;
        We only want to dilate values 1 and 3
        """

        # flush raster data to the disk
        slantrange_layover_shadow_mask_raster.close_dataset()
        del slantrange_layover_shadow_mask_raster

        # read layover/shadow mask
        gdal_ds = gdal.Open(path_layover_shadow_mask_file, gdal.GA_Update)
        gdal_band = gdal_ds.GetRasterBand(1)
        slantrange_layover_shadow_mask = gdal_band.ReadAsArray()

        # save layover pixels and substitute them with 0
        ind = np.where(slantrange_layover_shadow_mask == 2)
        slantrange_layover_shadow_mask[ind] = 0

        # perform grey dilation
        slantrange_layover_shadow_mask = ndimage.grey_dilation(
            slantrange_layover_shadow_mask, size=(shadow_dilation_size, shadow_dilation_size)
        )

        # restore layover pixels
        slantrange_layover_shadow_mask[ind] = 2

        # write dilated layover/shadow mask
        gdal_band.WriteArray(slantrange_layover_shadow_mask)

        # flush updates to the disk
        gdal_band.FlushCache()
        gdal_band = None
        gdal_ds = None

        slantrange_layover_shadow_mask_raster = isce3.io.Raster(path_layover_shadow_mask_file)

    # geocode the layover/shadow mask
    geo = isce3.geocode.GeocodeFloat32()
    geo.orbit = orbit
    geo.ellipsoid = ellipsoid
    geo.doppler = doppler
    geo.threshold_geo2rdr = threshold_geo2rdr
    geo.numiter_geo2rdr = numiter_geo2rdr
    geo.data_interpolator = 'NEAREST'
    geo.geogrid(
        float(geogrid_in.start_x),
        float(geogrid_in.start_y),
        float(geogrid_in.spacing_x),
        float(geogrid_in.spacing_y),
        int(geogrid_in.width),
        int(geogrid_in.length),
        int(geogrid_in.epsg),
    )

    geocoded_layover_shadow_mask_raster = isce3.io.Raster(
        filename_out, geogrid_in.width, geogrid_in.length, 1, gdal.GDT_Byte, output_raster_format
    )

    if geocode_options is None:
        geocode_options = {}

    if memory_mode is not None:
        geocode_options['memory_mode'] = memory_mode

    geo.geocode(
        radar_grid=radar_grid,
        input_raster=slantrange_layover_shadow_mask_raster,
        output_raster=geocoded_layover_shadow_mask_raster,
        dem_raster=dem_raster,
        output_mode=isce3.geocode.GeocodeOutputMode.INTERP,
        **geocode_options,
    )

    # flush data to the disk
    geocoded_layover_shadow_mask_raster.close_dataset()
    del geocoded_layover_shadow_mask_raster

    return slantrange_layover_shadow_mask_raster


def _create_raster_obj(
    output_dir,
    product_id,
    layer_name,
    dtype,
    shape,
    radar_grid_file_dict,
    output_obj_list,
):
    """Create an ISCE3 raster object (GTiff) for a radar geometry layer.

    Parameters
    ----------
    output_dir: str
           Output directory
    product_id: str
           Product ID
    dtype:: gdal.DataType
           GDAL data type
    shape: list
           Shape of the output raster
    radar_grid_file_dict: dict
           Dictionary that will hold the name of the output file
           referenced by the contents of `ds_hdf5` (dict key)
    output_obj_list: list
           Mutable list of output raster objects

    Returns
    -------
    raster_obj : isce3.io.Raster
           ISCE3 raster object
    """
    ds_name = f'{product_id}_{layer_name}'
    output_file = os.path.join(output_dir, ds_name) + '.tif'
    raster_obj = isce3.io.Raster(output_file, shape[2], shape[1], shape[0], dtype, 'GTiff')
    output_obj_list.append(raster_obj)
    radar_grid_file_dict[layer_name] = output_file
    return raster_obj


def save_intermediate_geocode_files(
    geogrid,
    dem_interp_method_enum,
    product_id,
    output_dir,
    extension,
    dem_raster,
    radar_grid_file_dict,
    lookside,
    wavelength,
    orbit,
    doppler=None,
):
    if doppler is None:
        doppler = isce3.core.LUT2d()

    # FIXME: Computation of range slope is not merged to ISCE yet
    output_obj_list = []
    layers_nbands = 1
    shape = [layers_nbands, geogrid.length, geogrid.width]
    names = [
        LAYER_NAME_LOCAL_INCIDENCE_ANGLE,
        LAYER_NAME_INCIDENCE_ANGLE,
        LAYER_NAME_PROJECTION_ANGLE,
        LAYER_NAME_RTC_ANF_PROJECTION_ANGLE,
        # LAYER_NAME_RANGE_SLOPE, # FIXME
        LAYER_NAME_DEM,
    ]
    raster_objs = []
    for name in names:
        raster_obj = _create_raster_obj(
            output_dir,
            product_id,
            name,
            gdal.GDT_Float32,
            shape,
            radar_grid_file_dict,
            output_obj_list,
        )
        raster_objs.append(raster_obj)
    (
        local_incidence_angle_raster,
        incidence_angle_raster,
        projection_angle_raster,
        rtc_anf_projection_angle_raster,
        # range_slope_raster, # FIXME
        interpolated_dem_raster,
    ) = raster_objs

    # TODO review this (Doppler)!!!
    # native_doppler = burst.doppler.lut2d
    native_doppler = doppler
    native_doppler.bounds_error = False
    grid_doppler = doppler
    grid_doppler.bounds_error = False

    isce3.geogrid.get_radar_grid(
        lookside,
        wavelength,
        dem_raster,
        geogrid,
        orbit,
        native_doppler,
        grid_doppler,
        dem_interp_method_enum,
        incidence_angle_raster=incidence_angle_raster,
        local_incidence_angle_raster=local_incidence_angle_raster,
        projection_angle_raster=projection_angle_raster,
        simulated_radar_brightness_raster=rtc_anf_projection_angle_raster,
        interpolated_dem_raster=interpolated_dem_raster,
        # range_slope_angle_raster=range_slope_raster, # FIXME
    )
    for obj in output_obj_list:
        del obj


def rtc(slc, geogrid, opts):
    # Common initializations
    t_start = time.time()
    output_dir = str(opts.output_dir)
    product_id = slc.id
    os.makedirs(output_dir, exist_ok=True)

    raster_format = 'GTiff'
    raster_extension = 'tif'

    # Filenames
    geo_filename = f'{output_dir}/{product_id}.{raster_extension}'
    nlooks_file = f'{output_dir}/{product_id}_{LAYER_NAME_NUMBER_OF_LOOKS}.{raster_extension}'
    rtc_anf_file = f'{output_dir}/{product_id}_{opts.layer_name_rtc_anf}.{raster_extension}'
    rtc_anf_gamma0_to_sigma0_file = (
        f'{output_dir}/{product_id}_{LAYER_NAME_RTC_ANF_GAMMA0_TO_SIGMA0}.{raster_extension}'
    )
    radar_grid = slc.radar_grid
    orbit = slc.orbit
    wavelength = slc.wavelength
    lookside = radar_grid.lookside

    dem_raster = isce3.io.Raster(opts.dem_path)
    ellipsoid = isce3.core.Ellipsoid()
    doppler = slc.doppler_centroid_grid
    exponent = 2

    x_snap = geogrid.spacing_x
    y_snap = geogrid.spacing_y
    geogrid.start_x = np.floor(float(geogrid.start_x) / x_snap) * x_snap
    geogrid.start_y = np.ceil(float(geogrid.start_y) / y_snap) * y_snap

    # geocoding optional arguments
    geocode_kwargs = {}
    layover_shadow_mask_geocode_kwargs = {}

    if isinstance(slc, SicdSlc):
        input_filename = slc.filepath.parent / (slc.filepath.stem + '_beta0.tif')
        if not input_filename.exists():
            slc.create_complex_beta0(input_filename)
        input_filename = str(input_filename)
    elif isinstance(slc, S1BurstSlc):
        input_filename = slc.filepath.parent / (slc.filepath.stem + '_beta0.tif')
        if not input_filename.exists():
            slc.create_complex_beta0(input_filename, flag_thermal_correction=opts.apply_thermal_noise)
        input_filename = str(input_filename)
        sub_swaths = slc.apply_valid_data_masking()
        geocode_kwargs['sub_swaths'] = sub_swaths
        layover_shadow_mask_geocode_kwargs['sub_swaths'] = sub_swaths
    else:
        input_filename = str(slc.filepath)

    layover_shadow_mask_file = f'{output_dir}/{product_id}_{LAYER_NAME_LAYOVER_SHADOW_MASK}.{raster_extension}'
    logger.info(f'    computing layover shadow mask for {product_id}')
    radar_grid_layover_shadow_mask = radar_grid
    slantrange_layover_shadow_mask_raster = compute_layover_shadow_mask(
        radar_grid_layover_shadow_mask,
        orbit,
        geogrid,
        dem_raster,
        layover_shadow_mask_file,
        raster_format,
        output_dir,
        shadow_dilation_size=opts.shadow_dilation_size,
        threshold_rdr2geo=opts.rdr2geo_threshold,
        numiter_rdr2geo=opts.rdr2geo_numiter,
        threshold_geo2rdr=opts.geo2rdr_threshold,
        numiter_geo2rdr=opts.geo2rdr_numiter,
        memory_mode=opts.memory_mode_isce3,
        geocode_options=layover_shadow_mask_geocode_kwargs,
        doppler=doppler,
    )
    logger.info(f'file saved: {layover_shadow_mask_file}')
    if opts.apply_shadow_masking:
        geocode_kwargs['input_layover_shadow_mask_raster'] = slantrange_layover_shadow_mask_raster

    out_geo_nlooks_obj = isce3.io.Raster(nlooks_file, geogrid.width, geogrid.length, 1, gdal.GDT_Float32, raster_format)
    out_geo_rtc_obj = isce3.io.Raster(rtc_anf_file, geogrid.width, geogrid.length, 1, gdal.GDT_Float32, raster_format)
    out_geo_rtc_gamma0_to_sigma0_obj = isce3.io.Raster(
        rtc_anf_gamma0_to_sigma0_file, geogrid.width, geogrid.length, 1, gdal.GDT_Float32, raster_format
    )
    geocode_kwargs['out_geo_rtc_gamma0_to_sigma0'] = out_geo_rtc_gamma0_to_sigma0_obj
    if opts.apply_bistatic_delay or opts.apply_static_tropo:
        rg_lut, az_lut = compute_correction_lut(
            slc.source,
            dem_raster,
            output_dir,
            opts.correction_lut_range_spacing_in_meters,
            opts.correction_lut_azimuth_spacing_in_meters,
            opts.apply_bistatic_delay,
            opts.apply_static_tropo,
        )
        geocode_kwargs['az_time_correction'] = az_lut
        if rg_lut is not None:
            geocode_kwargs['slant_range_correction'] = rg_lut

    rdr_raster = isce3.io.Raster(input_filename)
    # Generate output geocoded burst raster
    geo_raster = isce3.io.Raster(
        geo_filename, geogrid.width, geogrid.length, rdr_raster.num_bands, gdal.GDT_Float32, raster_format
    )

    # init Geocode object depending on raster type
    if rdr_raster.datatype() == gdal.GDT_Float32:
        geo_obj = isce3.geocode.GeocodeFloat32()
    elif rdr_raster.datatype() == gdal.GDT_Float64:
        geo_obj = isce3.geocode.GeocodeFloat64()
    elif rdr_raster.datatype() == gdal.GDT_CFloat32:
        geo_obj = isce3.geocode.GeocodeCFloat32()
    elif rdr_raster.datatype() == gdal.GDT_CFloat64:
        geo_obj = isce3.geocode.GeocodeCFloat64()
    else:
        err_str = 'Unsupported raster type for geocoding'
        raise NotImplementedError(err_str)

    # init geocode members
    geo_obj.orbit = orbit
    geo_obj.ellipsoid = ellipsoid
    geo_obj.doppler = doppler
    geo_obj.threshold_geo2rdr = opts.geo2rdr_threshold
    geo_obj.numiter_geo2rdr = opts.geo2rdr_numiter

    # set data interpolator based on the geocode algorithm
    if opts.geocode_algorithm_isce3 == isce3.geocode.GeocodeOutputMode.INTERP:
        geo_obj.data_interpolator = opts.geocode_algorithm_isce3

    geo_obj.geogrid(
        geogrid.start_x,
        geogrid.start_y,
        geogrid.spacing_x,
        geogrid.spacing_y,
        geogrid.width,
        geogrid.length,
        geogrid.epsg,
    )

    geo_obj.geocode(
        radar_grid=radar_grid,
        input_raster=rdr_raster,
        output_raster=geo_raster,
        dem_raster=dem_raster,
        output_mode=opts.geocode_algorithm_isce3,
        geogrid_upsampling=opts.geogrid_upsampling,
        flag_apply_rtc=opts.apply_rtc,
        input_terrain_radiometry=opts.input_terrain_radiometry_isce3,
        output_terrain_radiometry=opts.terrain_radiometry_isce3,
        exponent=exponent,
        rtc_min_value_db=opts.rtc_min_value_db,
        rtc_upsampling=opts.rtc_upsampling,
        rtc_algorithm=opts.rtc_algorithm_isce3,
        abs_cal_factor=opts.abs_cal_factor,
        flag_upsample_radar_grid=opts.upsample_radar_grid,
        clip_min=opts.clip_min,
        clip_max=opts.clip_max,
        out_geo_nlooks=out_geo_nlooks_obj,
        out_geo_rtc=out_geo_rtc_obj,
        rtc_area_beta_mode=opts.rtc_area_beta_mode_isce3,
        # out_geo_rtc_gamma0_to_sigma0=out_geo_rtc_gamma0_to_sigma0_obj,
        input_rtc=None,
        output_rtc=None,
        dem_interp_method=opts.dem_interpolation_method_isce3,
        memory_mode=opts.memory_mode_isce3,
        **geocode_kwargs,
    )

    del geo_raster

    out_geo_nlooks_obj.close_dataset()
    del out_geo_nlooks_obj

    out_geo_rtc_obj.close_dataset()
    del out_geo_rtc_obj

    out_geo_rtc_gamma0_to_sigma0_obj.close_dataset()
    del out_geo_rtc_gamma0_to_sigma0_obj

    radar_grid_file_dict = {}
    save_intermediate_geocode_files(
        geogrid,
        opts.dem_interpolation_method_isce3,
        product_id,
        output_dir,
        raster_extension,
        dem_raster,
        radar_grid_file_dict,
        lookside,
        wavelength,
        orbit,
        doppler=doppler,
    )
    t_end = time.time()
    logger.info(f'elapsed time: {t_end - t_start}')


def pfa_prototype_geocode(sicd, geogrid, dem_path, output_dir):
    interp_method = isce3.core.DataInterpMethod.BIQUINTIC
    sigma0_data = sicd.load_scaled_data('sigma0', power=True)
    slc_lut = isce3.core.LUT2d(
        np.arange(sigma0_data.shape[1]), np.arange(sigma0_data.shape[0]), sigma0_data, interp_method
    )
    local2ecef = pyproj.Transformer.from_crs(f'EPSG:{geogrid.epsg}', 'EPSG:4978', always_xy=True)
    local2ll = pyproj.Transformer.from_crs(f'EPSG:{geogrid.epsg}', 'EPSG:4326', always_xy=True)
    dem_raster = isce3.io.Raster(str(dem_path))
    dem = isce3.geometry.DEMInterpolator()
    dem.load_dem(dem_raster)
    dem.interp_method = interp_method
    output = np.zeros((geogrid.length, geogrid.width), dtype=np.float32)
    mask = np.zeros((geogrid.length, geogrid.width), dtype=bool)

    n_iters = geogrid.width * geogrid.length
    for i, j in tqdm(itertools.product(range(geogrid.width), range(geogrid.length)), total=n_iters):
        x = geogrid.start_x + (i * geogrid.spacing_x)
        y = geogrid.start_y + (j * geogrid.spacing_y)
        lon, lat = local2ll.transform(x, y)
        hae = dem.interpolate_lonlat(np.deg2rad(x), np.deg2rad(y))  # ISCE3 expects lat/lon to be in radians!
        ecef_x, ecef_y, ecef_z = local2ecef.transform(x, y, hae)
        row, col = sicd.geo2rowcol(np.array([(ecef_x, ecef_y, ecef_z)]))[0]
        if slc_lut.contains(row, col):
            output[j, i] = slc_lut.eval(row, col)
            mask[j, i] = 1

    output[mask == 0] = np.nan
    output_path = output_dir / f'{sicd.id}_{sicd.polarization}.tif'
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(str(output_path), geogrid.width, geogrid.length, 1, gdal.GDT_Float32)
    # account for pixel as area
    start_x = geogrid.start_x - (geogrid.spacing_x / 2)
    start_y = geogrid.start_y + (geogrid.spacing_y / 2)
    out_ds.SetGeoTransform([start_x, geogrid.spacing_x, 0, start_y, 0, geogrid.spacing_y])
    out_ds.SetProjection(pyproj.CRS(geogrid.epsg).to_wkt())
    out_ds.GetRasterBand(1).WriteArray(output)
    out_ds.GetRasterBand(1).SetNoDataValue(np.nan)
    out_ds.SetMetadata({'AREA_OR_POINT': 'Area'})
    out_ds = None
