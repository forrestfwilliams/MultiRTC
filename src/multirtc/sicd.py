from datetime import timedelta
from pathlib import Path
from typing import Optional

import isce3
import numpy as np
import pyproj
from numpy.polynomial.polynomial import polyval2d
from osgeo import gdal
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Point, Polygon

from multirtc import define_geogrid
from multirtc.base import Slc, print_wkt, to_isce_datetime


def check_poly_order(poly):
    assert len(poly.Coefs) == poly.order1 + 1, 'Polynomial order does not match number of coefficients'


class SicdSlc:
    """Base class for SICD SLCs."""

    def __init__(self, sicd_path: Path):
        self.reader = SICDReader(str(sicd_path.expanduser().resolve()))
        sicd = self.reader.get_sicds_as_tuple()[0]
        self.source = sicd
        self.id = Path(sicd_path).with_suffix('').name
        self.filepath = Path(sicd_path)
        self.footprint = Polygon([(ic.Lon, ic.Lat) for ic in sicd.GeoData.ImageCorners])
        self.center = Point(sicd.GeoData.SCP.LLH.Lon, sicd.GeoData.SCP.LLH.Lat)
        self.local_epsg = define_geogrid.get_point_epsg(self.center.y, self.center.x)
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
        self.raw_time_coa_poly = sicd.Grid.TimeCOAPoly
        self.arp_pos = sicd.SCPCOA.ARPPos.get_array()
        self.scp_pos = sicd.GeoData.SCP.ECF.get_array()
        self.look_angle = int(sicd.SCPCOA.AzimAng + 180) % 360
        self.beta0 = sicd.Radiometric.BetaZeroSFPoly
        self.sigma0 = sicd.Radiometric.SigmaZeroSFPoly
        self.supports_bistatic_delay = False
        self.supports_static_tropo = False

    def get_xrow_ycol(
        self, rowrange: Optional[tuple] = None, colrange: Optional[tuple] = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate xrow and ycol index arrrays.

        Args:
            rowrange: Optional tuple specifying the range of rows (start, end).
            colrange: Optional tuple specifying the range of columns (start, end).

        Returns:
            Two 2D numpy arrays, xrow and ycol, representing the row and column indices
            adjusted by the SCP index and scaled by the spacing.
        """
        rowlen = self.shape[0] if rowrange is None else rowrange[1] - rowrange[0]
        collen = self.shape[1] if colrange is None else colrange[1] - colrange[0]
        rowoffset = self.scp_index[0] if rowrange is None else self.scp_index[0] + rowrange[0]
        coloffset = self.scp_index[1] if colrange is None else self.scp_index[1] + colrange[0]

        irow = np.tile(np.arange(rowlen), (collen, 1)).T
        irow -= rowoffset
        xrow = irow * self.spacing[0]

        icol = np.tile(np.arange(collen), (rowlen, 1))
        icol -= coloffset
        ycol = icol * self.spacing[1]
        return xrow, ycol

    def load_scaled_data(
        self, scale: str, power: bool = False, rowrange: Optional[tuple] = None, colrange: Optional[tuple] = None
    ) -> np.ndarray:
        """Load scaled data from the SICD file.

        Args:
            scale: Scale type, either 'beta0' or 'sigma0'.
            power: If True, return power (squared magnitude), otherwise return complex data.
            rowrange: Optional tuple specifying the range of rows (start, end).
            colrange: Optional tuple specifying the range of columns (start, end).

        Returns:
            2D numpy array of scaled data, either power or complex, based on the scale type.
        """
        if scale == 'beta0':
            coeff = self.beta0.Coefs
        elif scale == 'sigma0':
            coeff = self.sigma0.Coefs
        else:
            raise ValueError(f'Scale must be either "beta0" or "sigma0", got {scale}')

        xrow, ycol = self.get_xrow_ycol(rowrange=rowrange, colrange=colrange)
        if colrange is not None and rowrange is not None:
            data = self.reader[rowrange[0] : rowrange[1], colrange[0] : colrange[1]]
        elif colrange is None and rowrange is None:
            data = self.reader[:, :]
        else:
            raise ValueError('Both xrange and yrange must be provided or neither.')

        scale_factor = polyval2d(xrow, ycol, coeff)
        del xrow, ycol  # deleting for memory management

        if power:
            data = (data.real**2 + data.imag**2) * scale_factor
        else:
            data = data * np.sqrt(scale_factor)
        return data

    def create_complex_beta0(self, outpath: str, row_iter: int = 256) -> None:
        """Create a complex beta0 image from the SICD data.
        Calculates the beta0 data in chunks to avoid memory issues.

        Args:
            outpath: Path to save the output beta0 TIFF file.
            row_iter: Number of rows to process in each chunk.
        """
        driver = gdal.GetDriverByName('GTiff')
        # Shape transposed for ISCE3 expectations
        ds = driver.Create(str(outpath), self.shape[0], self.shape[1], 1, gdal.GDT_CFloat32)
        band = ds.GetRasterBand(1)
        n_chunks = int(np.floor(self.shape[0] // row_iter)) + 1
        for i in range(n_chunks):
            start_row = i * row_iter
            end_row = min((i + 1) * row_iter, self.shape[0])
            rowrange = [start_row, end_row]
            colrange = [0, self.shape[1]]
            scaled_data = self.load_scaled_data('beta0', power=False, rowrange=rowrange, colrange=colrange)
            # Shape transposed for ISCE3 expectations
            if self.az_reversed:
                scaled_data = scaled_data[:, ::-1].T
            else:
                scaled_data = scaled_data.T
            # Offset transposed to match ISCE3 expectations
            band.WriteArray(scaled_data, xoff=start_row, yoff=0)

        band.FlushCache()
        ds.FlushCache()
        ds = None


class SicdRzdSlc(Slc, SicdSlc):
    """Class for SICD SLCs with range zero doppler grids."""

    def __init__(self, sicd_path: Path):
        super().__init__(sicd_path)
        assert self.source.Grid.Type == 'RGZERO', 'Only range zero doppler grids are supported for by this class'
        first_col_time = self.source.RMA.INCA.TimeCAPoly(-self.shift[1] * self.spacing[1])
        last_col_time = self.source.RMA.INCA.TimeCAPoly((self.shape[1] - self.shift[1]) * self.spacing[1])
        self.az_reversed = last_col_time < first_col_time
        self.sensing_start = min(first_col_time, last_col_time)
        self.sensing_end = max(first_col_time, last_col_time)
        self.starting_range = self.get_starting_range(0)
        self.az_reversed = last_col_time < first_col_time
        self.prf = self.shape[1] / (self.sensing_end - self.sensing_start)
        self.orbit = self.get_orbit()
        self.radar_grid = self.get_radar_grid()
        self.doppler_centroid_grid = isce3.core.LUT2d()
        self.supports_rtc = True

    def get_starting_range(self, col: int) -> float:
        assert 0 <= col < self.shape[1], 'Row index out of bounds'
        ycol = (col - self.shift[1]) * self.spacing[1]
        xrow = -self.shift[0] * self.spacing[0]  # fixing to first row
        inca_time = self.source.RMA.INCA.TimeCAPoly(ycol)
        arp_pos = self.arp_pos_poly(inca_time)
        row_offset = self.source.Grid.Row.UVectECF.get_array() * xrow
        col_offset = self.source.Grid.Col.UVectECF.get_array() * ycol
        grid_pos = self.source.GeoData.SCP.ECF.get_array() + row_offset + col_offset
        starting_range = np.linalg.norm(arp_pos - grid_pos)
        return starting_range

    def get_orbit(self) -> isce3.core.Orbit:
        """Define the orbit for the SLC.

        Returns:
            An instance of isce3.core.Orbit representing the orbit.
        """
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

    def get_radar_grid(self) -> isce3.product.RadarGridParameters:
        """Define the radar grid parameters for the SLC.

        Returns:
            An instance of isce3.product.RadarGridParameters representing the radar grid.
        """
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

    def create_geogrid(self, spacing_meters: int) -> isce3.product.GeoGridParameters:
        return define_geogrid.generate_geogrids(self, spacing_meters, self.local_epsg)

    def _print_wkt(self):
        return print_wkt(self)


class SicdPfaSlc(Slc, SicdSlc):
    """Class for SICD SLCs with PFA (Polar Format Algorithm) grids."""

    def __init__(self, sicd_path: Path):
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
        # TOOD: this may not always be true, will need to figure out a way to check
        self.az_reversed = False
        # Without ISCE3 support for PFA grids, these properties are undefined
        self.starting_range = np.nan
        self.radar_grid = None
        self.doppler_centroid_grid = None
        self.prf = np.nan
        self.az_reversed = False
        self.supports_rtc = False

    def get_orbit(self) -> isce3.core.Orbit:
        """Define the orbit for the SLC.
        PFA data has a constant COA time, so we create a simple orbit

        Returns:
            An instance of isce3.core.Orbit representing the orbit.
        """
        svs = []
        sensing_start_isce = to_isce_datetime(self.scp_time)
        for offset_sec in range(-10, 10):
            t = self.scp_time + timedelta(offset_sec)
            t_isce = to_isce_datetime(t)
            pos = self.arp_vel * offset_sec + self.arp_pos
            svs.append(isce3.core.StateVector(t_isce, pos, self.arp_vel))
        return isce3.core.Orbit(svs, sensing_start_isce)

    def calculate_range_range_rate_offset(self) -> np.ndarray:
        """Calculate the range and range rate offset for PFA data.

        Returns:
            A 2D numpy array containing the range and range rate offsets.
        """
        arp_minus_scp = self.arp_pos - self.scp_pos
        range_scp_to_coa = np.linalg.norm(arp_minus_scp, axis=-1)
        range_rate_scp_to_coa = np.sum(self.arp_vel * arp_minus_scp, axis=-1) / range_scp_to_coa
        rrdot_offset = np.array([range_scp_to_coa, range_rate_scp_to_coa])
        return rrdot_offset

    def calculate_transform_matrix(self) -> np.ndarray:
        """Define the matrix for transforming PFA grid coordinates to range and range rate.

        Returns:
            A 2x2 numpy array representing the transformation matrix.
        """
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

    def create_geogrid(self, spacing_meters: int) -> isce3.product.GeoGridParameters:
        """Create a geogrid for the PFA SLC.
        Note: Unlike other Slc subclasses, the PFA geogrid is always defined in EPSG 4326 (Lat/Lon).

        Args:
            spacing_meters: Spacing in meters for the geogrid.

        Returns:
            isce3.product.GeoGridParameters: The generated geogrid parameters.
        """
        ecef = pyproj.CRS(4978)  # ECEF on WGS84 Ellipsoid
        lla = pyproj.CRS(4979)  # WGS84 lat/lon/ellipsoid height
        local_utm = pyproj.CRS(define_geogrid.get_point_epsg(self.center.y, self.center.x))
        lla2utm = pyproj.Transformer.from_crs(lla, local_utm, always_xy=True)
        utm2lla = pyproj.Transformer.from_crs(local_utm, lla, always_xy=True)
        ecef2lla = pyproj.Transformer.from_crs(ecef, lla, always_xy=True)

        lla_point = (self.center.x, self.center.y)
        utm_point = lla2utm.transform(*lla_point)
        utm_point_shift = (utm_point[0] + spacing_meters, utm_point[1])
        lla_point_shift = utm2lla.transform(*utm_point_shift)
        x_spacing = lla_point_shift[0] - lla_point[0]
        y_spacing = -1 * x_spacing

        points = np.array([(0, 0), (0, self.shape[1]), self.shape, (self.shape[0], 0)])
        geos = self.rowcol2geo(points, self.scp_hae)

        points = np.vstack(ecef2lla.transform(geos[:, 0], geos[:, 1], geos[:, 2])).T
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
            epsg=4326,
        )
        geogrid_snapped = define_geogrid.snap_geogrid(geogrid, geogrid.spacing_x, geogrid.spacing_y)
        return geogrid_snapped
