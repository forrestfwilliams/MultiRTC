from datetime import timedelta
from pathlib import Path

import isce3
import numpy as np
import pyproj
from numpy.polynomial.polynomial import polyval2d
from osgeo import gdal
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Point, Polygon

from multirtc import define_geogrid
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
        self.scp_hae = sicd.GeoData.SCP.LLH.HAE
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
        self.raw_time_coa_poly = sicd.Grid.TimeCOAPoly
        last_line_time = self.raw_time_coa_poly(0, self.shape[1] - self.shift[1])
        first_line_time = self.raw_time_coa_poly(0, -self.shift[1])
        self.az_reversed = last_line_time > first_line_time
        self.arp_pos = sicd.SCPCOA.ARPPos.get_array()
        self.scp_pos = sicd.GeoData.SCP.ECF.get_array()
        azimuth_angle, elevation_angle = self.calculate_look_angles()
        assert np.isclose(azimuth_angle, sicd.SCPCOA.AzimAng), 'Azimuth angle does not match SCPCOA AzimuthAng'
        assert np.isclose(elevation_angle, sicd.SCPCOA.GrazeAng), 'Elevation angle does not match SCPCOA ElevationAng'
        self.look_angle = int(azimuth_angle + 180) % 360
        self.beta0 = sicd.Radiometric.BetaZeroSFPoly
        self.sigma0 = sicd.Radiometric.SigmaZeroSFPoly

    def calculate_look_angles(self):
        # Convert observer ECEF to geodetic lat/lon/alt
        transformer = pyproj.Transformer.from_crs('EPSG:4978', 'EPSG:4979', always_xy=True)
        lon_deg, lat_deg, _ = transformer.transform(*self.scp_pos)

        # ECEF to ENU rotation matrix centered at the observer
        sin_lat = np.sin(np.deg2rad(lat_deg))
        cos_lat = np.cos(np.deg2rad(lat_deg))
        sin_lon = np.sin(np.deg2rad(lon_deg))
        cos_lon = np.cos(np.deg2rad(lon_deg))
        rotation_matrix = np.array(
            [
                [-sin_lon, cos_lon, 0],
                [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
                [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
            ]
        )

        # Transform the LOS vector into ENU frame
        los = np.subtract(self.arp_pos, self.scp_pos)
        topocentric = rotation_matrix @ los
        east, north, up = topocentric

        # Compute azimuth and elevation
        azimuth = np.arctan2(east, north) % (2 * np.pi)
        elevation = np.arcsin(up / np.linalg.norm(topocentric))
        return np.rad2deg(azimuth), np.rad2deg(elevation)

    def get_xrow_ycol(self) -> np.ndarray:
        """Calculate xrow and ycol SICD."""
        irow = np.tile(np.arange(self.shape[0]), (self.shape[1], 1)).T
        irow -= self.scp_index[0]
        xrow = irow * self.spacing[0]

        icol = np.tile(np.arange(self.shape[1]), (self.shape[0], 1))
        icol -= self.scp_index[1]
        ycol = icol * self.spacing[1]
        return xrow, ycol

    def load_data(self):
        reader = SICDReader(str(self.filepath))
        data = reader[:, :]
        return data

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
        assert self.source.Grid.Type == 'RGZERO', 'Only range zero doppler grids are supported for by this class'
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


class SicdPfaSlc(SlcTemplate, SicdSlc):
    def __init__(self, sicd_path):
        super().__init__(sicd_path)
        assert self.source.ImageFormation.ImageFormAlgo == 'PFA', 'Only PFA-focused data are supported by this class'
        assert self.source.Grid.Type == 'RGAZIM', 'Only range azimuth grids are supported by this class'
        assert self.raw_time_coa_poly.Coefs.size == 1, 'Only constant COA time is currently supported'
        self.coa_time = self.raw_time_coa_poly.Coefs[0][0]
        self.arp_vel = self.source.SCPCOA.ARPVel.get_array()
        self.scp_time = self.reference_time + timedelta(self.source.SCPCOA.SCPTime)
        self.sensing_start = self.coa_time
        self.pfa_vars = self.source.PFA
        self.orbit = self.get_orbit()
        self.rrdot_offset = self.calculate_range_range_rate_offset()
        self.transform_matrix = self.calculate_transform_matrix()
        self.transform_matrix_inv = np.linalg.inv(self.transform_matrix)
        # Without ISCE3 support for PFA grids, these properties are undefined
        self.radar_grid = None
        self.doppler_centroid_grid = None
        self.prf = np.nan

    def get_orbit(self):
        svs = []
        sensing_start_isce = to_isce_datetime(self.scp_time)
        for offset_sec in range(-10, 10):
            t = self.scp_time + timedelta(offset_sec)
            t_isce = to_isce_datetime(t)
            pos = self.arp_vel * offset_sec + self.arp_pos
            svs.append(isce3.core.StateVector(t_isce, pos, self.arp_vel))
        return isce3.core.Orbit(svs, sensing_start_isce)

    def calculate_range_range_rate_offset(self):
        arp_minus_scp = self.arp_pos - self.scp_pos
        range_scp_to_coa = np.linalg.norm(arp_minus_scp, axis=-1)
        range_rate_scp_to_coa = np.sum(self.arp_vel * arp_minus_scp, axis=-1) / range_scp_to_coa
        rrdot_offset = np.array([range_scp_to_coa, range_rate_scp_to_coa])
        return rrdot_offset

    def calculate_transform_matrix(self):
        polar_ang_poly = self.pfa_vars.PolarAngPoly
        spatial_freq_sf_poly = self.pfa_vars.SpatialFreqSFPoly
        polar_ang_poly_der = polar_ang_poly.derivative(der_order=1, return_poly=True)
        spatial_freq_sf_poly_der = spatial_freq_sf_poly.derivative(der_order=1, return_poly=True)

        polar_ang_poly_der = polar_ang_poly.derivative(der_order=1, return_poly=True)
        spatial_freq_sf_poly_der = spatial_freq_sf_poly.derivative(der_order=1, return_poly=True)

        thetaTgtCoa = polar_ang_poly(self.coa_time)
        dThetaDtTgtCoa = polar_ang_poly_der(self.coa_time)

        # Compute polar aperture scale factor (KSF) and derivative
        # wrt polar angle
        ksfTgtCoa = spatial_freq_sf_poly(thetaTgtCoa)
        dKsfDThetaTgtCoa = spatial_freq_sf_poly_der(thetaTgtCoa)

        # Compute spatial frequency domain phase slopes in Ka and Kc directions
        # NB: sign for the phase may be ignored as it is cancelled
        # in a subsequent computation.
        dPhiDKaTgtCoa = np.array([np.cos(thetaTgtCoa), np.sin(thetaTgtCoa)])
        dPhiDKcTgtCoa = np.array([-np.sin(thetaTgtCoa), np.cos(thetaTgtCoa)])

        transform_matrix = np.zeros((2, 2))
        transform_matrix[0, :] = ksfTgtCoa * dPhiDKaTgtCoa
        transform_matrix[1, :] = dThetaDtTgtCoa * (dKsfDThetaTgtCoa * dPhiDKaTgtCoa + ksfTgtCoa * dPhiDKcTgtCoa)
        return transform_matrix

    def rowcol2geo(self, rc: np.ndarray, hae: float) -> np.ndarray:
        """Transform grid (row, col) coordinates to ECEF coordinates.

        Args:
            rc: 2D array of (row, col) coordinates
            hae: Height above ellipsoid (meters)

        Returns:
            np.ndarray: ECEF coordinates
        """
        dem = isce3.geometry.DEMInterpolator(hae)
        elp = isce3.core.Ellipsoid()
        rgaz = (rc - np.array(self.shift)[None, :]) * np.array(self.spacing)[None, :]
        rrdot = np.dot(self.transform_matrix, rgaz.T) + self.rrdot_offset[:, None]
        side = isce3.core.LookSide(1) if self.lookside == 'left' else isce3.core.LookSide(-1)
        pts_ecf = []
        wvl = 1.0
        for pt in rrdot.T:
            r = pt[0]
            dop = -pt[1] * 2 / wvl
            llh = isce3.geometry.rdr2geo(0.0, r, self.orbit, side, dop, wvl, dem, threshold=1.0e-8, maxiter=50)
            pts_ecf.append(elp.lon_lat_to_xyz(llh))
        return np.vstack(pts_ecf)

    def geo2rowcol(self, xyz: np.ndarray) -> np.ndarray:
        """Transform ECEF xyz to (row, col).

        Args:
            xyz: ECEF coordinates

        Returns:
            (row, col) coordinates
        """
        rrdot = np.zeros((2, xyz.shape[0]))
        rrdot[0, :] = np.linalg.norm(xyz - self.arp_pos[None, :], axis=1)
        rrdot[1, :] = np.dot(-self.arp_vel, (xyz - self.arp_pos[None, :]).T) / rrdot[0, :]
        rgaz = np.dot(self.transform_matrix_inv, (rrdot - self.rrdot_offset[:, None]))
        rgaz /= np.array(self.spacing)[:, None]
        rgaz += np.array(self.shift)[:, None]
        row_col = rgaz.T.copy()
        return row_col

    def create_geogrid(self, spacing_meters):
        epsg_local = define_geogrid.get_point_epsg(self.center.y, self.center.x)
        ecef = pyproj.CRS(4978)  # ECEF on WGS84 Ellipsoid
        local = pyproj.CRS(epsg_local)
        ecef2local = pyproj.Transformer.from_crs(ecef, local, always_xy=True)
        x_spacing = spacing_meters
        y_spacing = -1 * spacing_meters

        points = np.array([(0, 0), (0, self.shape[1]), self.shape, (self.shape[0], 0)])
        geos = self.rowcol2geo(points, hae=self.scp_hae)

        points = np.vstack(ecef2local.transform(geos[:, 0], geos[:, 1], geos[:, 2])).T
        minx, maxx = np.min(points[:, 0]), np.max(points[:, 0])
        miny, maxy = np.min(points[:, 1]), np.max(points[:, 1])

        width = (maxx - minx) // x_spacing
        length = (maxy - miny) // np.abs(y_spacing)
        geogrid = isce3.product.GeoGridParameters(
            start_x=float(minx),
            start_y=float(maxy),
            spacing_x=float(x_spacing),
            spacing_y=float(y_spacing),
            length=int(length),
            width=int(width),
            epsg=epsg_local,
        )
        return geogrid
