from pathlib import Path

import isce3
import numpy as np
import s1reader
from osgeo import gdal

from multirtc import define_geogrid
from multirtc.base import Slc, from_isce_datetime, to_isce_datetime


gdal.UseExceptions()


class S1BurstSlc(Slc):
    """Class representing a Sentinel-1 burst SLC (Single Look Complex) product."""

    def __init__(self, safe_path: Path, orbit_path: Path, burst_name: str):
        _, burst_id, swath, _, polarization, _ = burst_name.split('_')
        burst_id = int(burst_id)
        swath_num = int(swath[2])
        bursts = s1reader.load_bursts(str(safe_path), str(orbit_path), swath_num, polarization)
        burst = [b for b in bursts if str(b.burst_id).endswith(f'{burst_id}_{swath.lower()}')][0]
        del bursts
        vrt_path = safe_path.parent / f'{burst_name}.vrt'
        burst.slc_to_vrt_file(vrt_path)
        self.id = burst_name
        self.filepath = vrt_path
        self.footprint = burst.border[0]
        self.center = burst.center
        self.local_epsg = define_geogrid.get_point_epsg(self.center.y, self.center.x)
        self.lookside = 'right'
        self.wavelength = burst.wavelength
        self.polarization = burst.polarization
        self.shape = burst.shape
        self.range_pixel_spacing = burst.range_pixel_spacing
        self.reference_time = from_isce_datetime(burst.orbit.reference_epoch)
        self.sensing_start = (burst.sensing_start - self.reference_time).total_seconds()
        self.starting_range = burst.starting_range
        self.prf = 1 / burst.azimuth_time_interval
        self.orbit = burst.orbit
        self.doppler_centroid_grid = isce3.core.LUT2d()
        self.radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start,
            wavelength=self.wavelength,
            prf=self.prf,
            starting_range=self.starting_range,
            range_pixel_spacing=self.range_pixel_spacing,
            lookside=isce3.core.LookSide.Right,
            length=self.shape[0],
            width=self.shape[1],
            ref_epoch=to_isce_datetime(self.reference_time),
        )
        self.first_valid_line = burst.first_valid_line
        self.last_valid_line = burst.last_valid_line
        self.first_valid_sample = burst.first_valid_sample
        self.last_valid_sample = burst.last_valid_sample
        self.source = burst
        self.supports_rtc = True
        self.supports_bistatic_delay = True
        self.supports_static_tropo = True

    def create_geogrid(self, spacing_meters: int) -> isce3.product.GeoGridParameters:
        return define_geogrid.generate_geogrids(self, spacing_meters, self.local_epsg)

    def apply_valid_data_masking(self) -> isce3.product.SubSwaths:
        """Extract burst boundaries and create sub_swaths object to mask invalid radar samples.

        Returns:
           SubSwaths object with valid samples set according to the burst boundaries.
        """
        n_subswaths = 1
        sub_swaths = isce3.product.SubSwaths(self.radar_grid.length, self.radar_grid.width, n_subswaths)
        last_range_sample = min([self.last_valid_sample, self.radar_grid.width])
        valid_samples_sub_swath = np.repeat(
            [[self.first_valid_sample, last_range_sample + 1]], self.radar_grid.length, axis=0
        )
        for i in range(self.first_valid_line):
            valid_samples_sub_swath[i, :] = 0
        for i in range(self.last_valid_line, self.radar_grid.length):
            valid_samples_sub_swath[i, :] = 0

        sub_swaths.set_valid_samples_array(1, valid_samples_sub_swath)
        return sub_swaths

    def create_complex_beta0(self, outpath: Path, flag_thermal_correction: bool = True) -> None:
        """Apply conversion to beta0 and optionally apply a thermal noise correction.

        Args:
            outpath: Path to save the corrected beta0 image.
            flag_thermal_correction: If True, apply thermal noise correction using the LUT from the source burst.
        """
        # Load the SLC of the burst
        slc_gdal_ds = gdal.Open(str(self.filepath))
        arr_slc_from = slc_gdal_ds.ReadAsArray()

        # Apply thermal noise correction
        if flag_thermal_correction:
            corrected_image = np.abs(arr_slc_from) ** 2 - self.source.thermal_noise_lut
            min_backscatter = 0
            max_backscatter = None
            corrected_image = np.clip(corrected_image, min_backscatter, max_backscatter)
        else:
            corrected_image = np.abs(arr_slc_from) ** 2

        # Apply absolute radiometric correction
        corrected_image = corrected_image / self.source.burst_calibration.beta_naught**2

        factor_mag = np.sqrt(corrected_image) / np.abs(arr_slc_from)
        factor_mag[np.isnan(factor_mag)] = 0.0
        corrected_image = arr_slc_from * factor_mag
        dtype = gdal.GDT_CFloat32

        # Save the corrected image
        drvout = gdal.GetDriverByName('GTiff')
        raster_out = drvout.Create(outpath, self.shape[1], self.shape[0], 1, dtype)
        band_out = raster_out.GetRasterBand(1)
        band_out.WriteArray(corrected_image)
        band_out.FlushCache()
        del band_out
