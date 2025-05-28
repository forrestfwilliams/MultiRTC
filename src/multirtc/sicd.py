from datetime import timedelta
from pathlib import Path

import isce3
import numpy as np
from numpy.polynomial.polynomial import polyval2d
from osgeo import gdal
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Point, Polygon

from multirtc.base import SlcTemplate, to_isce_datetime


def check_poly_order(poly):
    assert len(poly.Coefs) == poly.order1 + 1, 'Polynomial order does not match number of coefficients'


class SicdSlc:
    def __init__(self, sicd_path):
        reader = SICDReader(str(sicd_path.expanduser().resolve()))
        sicd = reader.get_sicds_as_tuple()[0]
        self.source = sicd
        self.id = Path(sicd_path).with_suffix('').name
        self.filepath = Path(sicd_path)
        self.footprint = Polygon([(ic.Lon, ic.Lat) for ic in sicd.GeoData.ImageCorners])
        self.center = Point(sicd.GeoData.SCP.LLH.Lon, sicd.GeoData.SCP.LLH.Lat)
        self.lookside = 'right' if sicd.SCPCOA.SideOfTrack == 'R' else 'left'

        center_frequency = sicd.RadarCollection.TxFrequency.Min + sicd.RadarCollection.TxFrequency.Max / 2
        self.wavelength = isce3.core.speed_of_light / center_frequency
        self.polarization = sicd.RadarCollection.RcvChannels[0].TxRcvPolarization.replace(':', '')
        self.shape = (sicd.ImageData.NumRows, sicd.ImageData.NumCols)
        self.spacing = (sicd.Grid.Row.SS, sicd.Grid.Col.SS)
        self.scp_index = (sicd.ImageData.SCPPixel.Row, sicd.ImageData.SCPPixel.Col)
        self.range_pixel_spacing = sicd.Grid.Row.SS
        self.reference_time = sicd.Timeline.CollectStart.item()
        self.shift = (
            sicd.ImageData.SCPPixel.Row - sicd.ImageData.FirstRow,
            sicd.ImageData.SCPPixel.Col - sicd.ImageData.FirstCol,
        )
        self.arp_pos_poly = sicd.Position.ARPPoly
        starting_row_pos = (
            sicd.GeoData.SCP.ECF.get_array()
            + sicd.Grid.Row.UVectECF.get_array() * (0 - self.shift[0]) * self.spacing[0]
        )
        self.starting_range = np.linalg.norm(sicd.SCPCOA.ARPPos.get_array() - starting_row_pos)
        last_line_time = sicd.Grid.TimeCOAPoly(0, self.shape[1] - self.shift[1])
        first_line_time = sicd.Grid.TimeCOAPoly(0, -self.shift[1])
        self.az_reversed = last_line_time >= first_line_time
        self.beta0 = sicd.Radiometric.BetaZeroSFPoly
        self.sigma0 = sicd.Radiometric.SigmaZeroSFPoly

    def get_xrow_ycol(self) -> np.ndarray:
        """Calculate xrow and ycol SICD."""
        irow = np.tile(np.arange(self.shape[0]), (self.shape[1], 1)).T
        irow -= self.scp_index[0]
        xrow = irow * self.row_mult

        icol = np.tile(np.arange(self.shape[1]), (self.shape[0], 1))
        icol -= self.scp_index[1]
        ycol = icol * self.col_mult
        return xrow, ycol

    def load_data(self):
        return self.source[:, :]

    def load_scaled_data(self, scale, power=False):
        if scale == 'beta0':
            coeff = self.beta0.Coefs
        elif scale == 'sigma0':
            coeff = self.sigma0.Coefs
        else:
            raise ValueError(f'Scale must be either "beta0" or "sigma0", got {scale}')

        xrow, ycol = self.get_xrow_ycol()
        scale_factor = polyval2d(xrow, ycol, coeff)
        data = self.load_data()
        if power:
            scaled_data = (data.real**2 + data.imag**2) * scale_factor
        else:
            scaled_data = data * np.sqrt(scale_factor)
        return scaled_data

    def write_complex_beta0(self, outpath, isce_format=True):
        scaled_data = self.load_scaled_data('beta0', power=False)
        if isce_format:
            if self.az_reversed:
                scaled_data = scaled_data[:, ::-1].T
            else:
                scaled_data = scaled_data.T

        driver = gdal.GetDriverByName('GTiff')
        ds = driver.Create(str(outpath), scaled_data.shape[1], scaled_data.shape[0], 1, gdal.GDT_CFloat32)
        band = ds.GetRasterBand(1)
        band.WriteArray(scaled_data)
        band.FlushCache()
        ds = None

    def create_complex_beta0(self, outpath, isce_format=True):
        xrow, ycol = self.get_xrow_ycol()
        scale_factor = np.sqrt(polyval2d(xrow, ycol, self.beta0_coeff))
        data = self.load_data()
        scaled_data = data * scale_factor

        if isce_format:
            if self.az_reversed:
                scaled_data = scaled_data[:, ::-1].T
            else:
                scaled_data = scaled_data.T

        driver = gdal.GetDriverByName('GTiff')
        ds = driver.Create(str(outpath), scaled_data.shape[1], scaled_data.shape[0], 1, gdal.GDT_CFloat32)
        band = ds.GetRasterBand(1)
        band.WriteArray(scaled_data)
        band.FlushCache()
        ds = None


class SicdRzdSlc(SlcTemplate, SicdSlc):
    def __init__(self, sicd_path):
        super().__init__(sicd_path)
        assert self.source.Grid.Type == 'RGZERO', 'Only range zero doppler grids supported for Capella data'
        first_col_time = self.source.RMA.INCA.TimeCAPoly(0 - self.shift[1])
        last_col_time = self.source.RMA.INCA.TimeCAPoly(self.shape[1] - self.shift[1])
        self.sensing_start = min(first_col_time, last_col_time)
        self.sensing_end = max(first_col_time, last_col_time)
        self.prf = self.shape[1] / (self.sensing_end - self.sensing_start)
        self.orbit = self.get_orbit()
        self.radar_grid = self.get_radar_grid()
        self.doppler_centroid_grid = isce3.core.LUT2d()

    def get_orbit(self):
        svs = []
        orbit_start = np.floor(self.sensing_start) - 10
        orbit_end = np.ceil(self.sensing_end) + 10
        for offset_sec in np.arange(orbit_start, orbit_end + 1, 1):
            t = self.sensing_start + offset_sec
            pos = self.arp_pos_poly(t)
            vel = self.arp_pos_poly.derivative_eval(t)
            t_isce = to_isce_datetime(self.reference_time + timedelta(seconds=t))
            svs.append(isce3.core.StateVector(t_isce, pos, vel))
        return isce3.core.Orbit(svs, to_isce_datetime(self.reference_time))

    def get_radar_grid(self):
        radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start,
            wavelength=self.wavelength,
            prf=self.prf,
            starting_range=self.starting_range,
            range_pixel_spacing=self.range_pixel_spacing,
            lookside=isce3.core.LookSide.Right if self.lookside == 'right' else isce3.core.LookSide.Left,
            length=self.shape[1],  # flipped for "shadows down" convention
            width=self.shape[0],  # flipped for "shadows down" convention
            ref_epoch=to_isce_datetime(self.reference_time),
        )
        return radar_grid
