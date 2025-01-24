from dataclasses import dataclass

import isce3
import numpy as np


@dataclass
class RtcOptions:
    output_dir: str
    dem_path: str
    apply_rtc: bool = True
    apply_thermal_noise: bool = True
    apply_abs_rad: bool = True
    apply_bistatic_delay: bool = True
    apply_static_tropo: bool = True
    apply_valid_samples_sub_swath_masking: bool = True
    apply_shadow_masking: bool = True
    dem_interpolation_method: str = 'biquintic'
    geocode_algorithm: str = 'area_projection'  # 'area_projection' or 'interp'
    correction_lut_azimuth_spacing_in_meters: int = 120
    correction_lut_range_spacing_in_meters: int = 120
    memory_mode: str = 'single_block'
    geogrid_upsampling: int = 1
    shadow_dilation_size: int = 0
    abs_cal_factor: int = 1
    clip_min: float = np.nan
    clip_max: float = np.nan
    upsample_radar_grid: bool = False
    terrain_radiometry: str = 'gamma0'  # 'gamma0' or 'sigma0'
    rtc_algorithm_type: str = 'area_projection'  # 'area_projection' or 'bilinear_distribution'
    input_terrain_radiometry: str = 'beta0'
    rtc_min_value_db: int = -30.0
    rtc_upsampling: int = 2
    rtc_area_beta_mode: str = 'auto'
    geo2rdr_threshold: float = 1.0e-7
    geo2rdr_numiter: int = 50
    rdr2geo_threshold: float = 1.0e-7
    rdr2geo_numiter: int = 25
    output_epsg: int = None
    x_spacing: int = 30
    y_spacing: int = 30

    def __post_init__(self):
        if not self.apply_rtc:
            if self.save_rtc_anf:
                raise ValueError('RTC ANF flags are only available with RTC enabled')
            if self.save_rtc_anf_gamma0_to_sigma0:
                raise ValueError('RTC ANF gamma0 to sigma0 flags are only available with RTC enabled')

        if self.terrain_radiometry == 'sigma0' and self.save_rtc_anf_gamma0_to_sigma0:
            raise ValueError('RTC ANF gamma0 to sigma0 flags are only available with output type set to gamma0')

        if self.apply_rtc:
            self.layer_name_rtc_anf = f'rtc_anf_{self.terrain_radiometry}_to_{self.input_terrain_radiometry}'
        else:
            self.layer_name_rtc_anf = ''

        if self.dem_interpolation_method == 'biquintic':
            self.dem_interpolation_method_isce3 = isce3.core.DataInterpMethod.BIQUINTIC
        else:
            raise ValueError(f'Invalid DEM interpolation method: {self.dem_interpolation_method}')

        if self.geocode_algorithm == 'area_projection':
            self.geocode_algorithm_isce3 = isce3.geocode.GeocodeOutputMode.AREA_PROJECTION
        elif self.geocode_algorithm == 'interp':
            self.geocode_algorithm_isce3 = isce3.geocode.GeocodeOutputMode.INTERP
        else:
            raise ValueError(f'Invalid geocode algorithm: {self.geocode_algorithm}')

        if self.memory_mode == 'single_block':
            self.memory_mode_isce3 = isce3.core.GeocodeMemoryMode.SingleBlock
        else:
            raise ValueError(f'Invalid memory mode: {self.memory_mode}')

        if self.rtc_algorithm_type == 'bilinear_distribution':
            self.rtc_algorithm_isce3 = isce3.geometry.RtcAlgorithm.RTC_BILINEAR_DISTRIBUTION
        elif self.rtc_algorithm_type == 'area_projection':
            self.rtc_algorithm_isce3 = isce3.geometry.RtcAlgorithm.RTC_AREA_PROJECTION
        else:
            raise ValueError(f'Invalid RTC algorithm: {self.rtc_algorithm_type}')

        if self.input_terrain_radiometry == 'sigma0':
            self.input_terrain_radiometry_isce3 = isce3.geometry.RtcInputTerrainRadiometry.SIGMA_NAUGHT_ELLIPSOID
        elif self.input_terrain_radiometry == 'beta0':
            self.input_terrain_radiometry_isce3 = isce3.geometry.RtcInputTerrainRadiometry.BETA_NAUGHT
        else:
            raise ValueError(f'Invalid input terrain radiometry: {self.input_terrain_radiometry}')

        if self.terrain_radiometry == 'sigma0':
            self.terrain_radiometry_isce3 = isce3.geometry.RtcOutputTerrainRadiometry.SIGMA_NAUGHT
        elif self.terrain_radiometry == 'gamma0':
            self.terrain_radiometry_isce3 = isce3.geometry.RtcOutputTerrainRadiometry.GAMMA_NAUGHT
        else:
            raise ValueError(f'Invalid terrain radiometry: {self.terrain_radiometry}')

        if self.apply_rtc:
            self.layer_name_rtc_anf = f'rtc_anf_{self.terrain_radiometry}_to_{self.input_terrain_radiometry}'
        else:
            self.layer_name_rtc_anf = ''

        if self.rtc_area_beta_mode == 'pixel_area':
            self.rtc_area_beta_mode_isce3 = isce3.geometry.RtcAreaBetaMode.PIXEL_AREA
        elif self.rtc_area_beta_mode == 'projection_angle':
            self.rtc_area_beta_mode_isce3 = isce3.geometry.RtcAreaBetaMode.PROJECTION_ANGLE
        elif self.rtc_area_beta_mode == 'auto' or self.rtc_area_beta_mode is None:
            self.rtc_area_beta_mode_isce3 = isce3.geometry.RtcAreaBetaMode.AUTO
        else:
            raise ValueError(f'Invalid area beta mode: {self.rtc_area_beta_mode}')

        self.x_snap = self.x_spacing
        self.y_snap = self.y_spacing
