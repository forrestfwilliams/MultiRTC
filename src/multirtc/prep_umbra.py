from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import isce3
import numpy as np
import pyproj
from numpy.polynomial.polynomial import polyval2d
from sarpy.io.complex.sicd import SICDReader
from shapely.geometry import Polygon

from multirtc import dem
from multirtc.define_geogrid import get_point_epsg


def check_poly_order(poly):
    assert len(poly.Coefs) == poly.order1 + 1, 'Polynomial order does not match number of coefficients'


@dataclass
class Point:
    x: float
    y: float


@dataclass
class UmbraSICD:
    id: str
    file_path: Path
    footprint: Polygon
    wavelength: float
    polarization: str
    lookside: str  # 'right' or 'left'
    shape: tuple[int, int]
    scp_index: tuple[int, int]
    scp_time: float
    scp_pos: np.ndarray
    center: Point
    scp_hae: float
    coa_time: float
    arp_pos: np.ndarray
    arp_vel: np.ndarray
    grid_shift: np.ndarray
    grid_mult: np.ndarray
    rrdot_offset: np.ndarray
    transform_matrix: np.ndarray
    transform_matrix_inv: np.ndarray
    orbit: isce3.core.Orbit
    beta0_coeff: np.ndarray
    sigma0_coeff: np.ndarray

    @staticmethod
    def calculate_orbit(sensing_start, pos_arp, vel_arp):
        svs = []
        sensing_start_isce = isce3.core.DateTime(datetime.utcfromtimestamp(sensing_start.astype(int) * 1e-9))
        for offset_sec in range(-5, 6):
            t = sensing_start + np.timedelta64(offset_sec, 's')
            t = isce3.core.DateTime(datetime.utcfromtimestamp(t.astype(int) * 1e-9))
            pos = vel_arp * offset_sec + pos_arp
            svs.append(isce3.core.StateVector(t, pos, vel_arp))
        return isce3.core.Orbit(svs, sensing_start_isce)

    @staticmethod
    def calculate_range_range_rate_offset(scp_pos, arp_pos, arp_vel, time_coa):
        arp_minus_scp = arp_pos - scp_pos
        range_scp_to_coa = np.linalg.norm(arp_minus_scp, axis=-1)
        range_rate_scp_to_coa = np.sum(arp_vel * arp_minus_scp, axis=-1) / range_scp_to_coa
        rrdot_offset = np.array([range_scp_to_coa, range_rate_scp_to_coa])
        return rrdot_offset

    @staticmethod
    def calculate_transform_matrix(pfa, time_coa):
        polar_ang_poly = pfa.PolarAngPoly
        spatial_freq_sf_poly = pfa.SpatialFreqSFPoly
        polar_ang_poly_der = polar_ang_poly.derivative(der_order=1, return_poly=True)
        spatial_freq_sf_poly_der = spatial_freq_sf_poly.derivative(der_order=1, return_poly=True)

        polar_ang_poly_der = polar_ang_poly.derivative(der_order=1, return_poly=True)
        spatial_freq_sf_poly_der = spatial_freq_sf_poly.derivative(der_order=1, return_poly=True)

        thetaTgtCoa = polar_ang_poly(time_coa)
        dThetaDtTgtCoa = polar_ang_poly_der(time_coa)

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

    @classmethod
    def from_sarpy_sicd(cls, sicd, file_path):
        center_frequency = sicd.RadarCollection.TxFrequency.Min + sicd.RadarCollection.TxFrequency.Max / 2
        wavelength = isce3.core.speed_of_light / center_frequency
        polarization = sicd.RadarCollection.RcvChannels[0].TxRcvPolarization.replace(':', '')
        lookside = 'right' if sicd.SCPCOA.SideOfTrack == 'R' else 'left'
        footprint = Polygon([(ic.Lon, ic.Lat) for ic in sicd.GeoData.ImageCorners])

        coef_time_coa = sicd.Grid.TimeCOAPoly.Coefs
        assert coef_time_coa.size == 1, 'Only constant COA time is currently supported'
        coa_time = coef_time_coa[0][0]

        arp_pos = np.array([sicd.SCPCOA.ARPPos.X, sicd.SCPCOA.ARPPos.Y, sicd.SCPCOA.ARPPos.Z])
        arp_vel = np.array([sicd.SCPCOA.ARPVel.X, sicd.SCPCOA.ARPVel.Y, sicd.SCPCOA.ARPVel.Z])
        scp_time = sicd.Timeline.CollectStart + np.timedelta64(int(sicd.SCPCOA.SCPTime * 1e9), 'ns')
        orbit = cls.calculate_orbit(scp_time, arp_pos, arp_vel)
        grid_shift = np.array(
            [
                sicd.ImageData.SCPPixel.Row - sicd.ImageData.FirstRow,
                sicd.ImageData.SCPPixel.Col - sicd.ImageData.FirstCol,
            ]
        )
        grid_mult = np.array([sicd.Grid.Row.SS, sicd.Grid.Col.SS])
        scp_pos = np.array([sicd.GeoData.SCP.ECF.X, sicd.GeoData.SCP.ECF.Y, sicd.GeoData.SCP.ECF.Z])
        scp_hae = sicd.GeoData.SCP.LLH.HAE
        rrdot_offset = cls.calculate_range_range_rate_offset(scp_pos, arp_pos, arp_vel, coa_time)
        transform_matrix = cls.calculate_transform_matrix(sicd.PFA, coa_time)
        beta0_coeff = sicd.Radiometric.BetaZeroSFPoly.Coefs
        sigma0_coeff = sicd.Radiometric.SigmaZeroSFPoly.Coefs
        umbra_sicd = cls(
            id=Path(file_path).with_suffix('').name,
            file_path=file_path,
            footprint=footprint,
            shape=(sicd.ImageData.NumRows, sicd.ImageData.NumCols),
            wavelength=wavelength,
            polarization=polarization,
            lookside=lookside,
            scp_index=(sicd.ImageData.SCPPixel.Row, sicd.ImageData.SCPPixel.Col),
            scp_time=scp_time,
            scp_pos=scp_pos,
            center=Point(sicd.GeoData.SCP.LLH.Lon, sicd.GeoData.SCP.LLH.Lat),
            scp_hae=scp_hae,
            coa_time=coa_time,
            arp_pos=arp_pos,
            arp_vel=arp_vel,
            grid_shift=grid_shift,
            grid_mult=grid_mult,
            rrdot_offset=rrdot_offset,
            transform_matrix=transform_matrix,
            transform_matrix_inv=np.linalg.inv(transform_matrix),
            orbit=orbit,
            beta0_coeff=beta0_coeff,
            sigma0_coeff=sigma0_coeff,
        )
        return umbra_sicd

    def rowcol2geo(self, rc: np.ndarray, hae: float = None) -> np.ndarray:
        """Transforma (row, col) coordinates to ECEF coordinates.

        Args:
            rc: Tuple of (row, col) coordinates
            hae: Height above ellipsoid (meters)

        Returns:
            np.ndarray: ECEF coordinates
        """
        if hae is None:
            hae = self.scp_hae

        dem = isce3.geometry.DEMInterpolator(hae)
        elp = isce3.core.Ellipsoid()
        rgaz = (rc - self.grid_shift[None, :]) * self.grid_mult[None, :]
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
        rgaz /= self.grid_mult[:, None]
        rgaz += self.grid_shift[:, None]
        row_col = rgaz.T.copy()
        return row_col

    def get_geogrid(self, x_spacing, y_spacing):
        ecef = pyproj.CRS(4978)  # ECEF on WGS84 Ellipsoid
        lla = pyproj.CRS(4979)  # WGS84 lat/lon/ellipsoid height
        local_utm = pyproj.CRS(get_point_epsg(self.center.y, self.center.x))
        lla2utm = pyproj.Transformer.from_crs(lla, local_utm, always_xy=True)
        utm2lla = pyproj.Transformer.from_crs(local_utm, lla, always_xy=True)
        ecef2lla = pyproj.Transformer.from_crs(ecef, lla, always_xy=True)

        lla_point = (self.center.x, self.center.y)
        utm_point = lla2utm.transform(*lla_point)
        utm_point_shift = (utm_point[0] + x_spacing, utm_point[1])
        lla_point_shift = utm2lla.transform(*utm_point_shift)
        x_spacing = lla_point_shift[0] - lla_point[0]
        y_spacing = -1 * x_spacing

        points = np.array([(0, 0), (0, self.shape[1]), self.shape, (self.shape[0], 0)])
        geos = self.rowcol2geo(points)

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
        return geogrid

    def get_xrow_ycol(self) -> np.ndarray:
        """Calculate xrow and ycol for the umbra_sicd."""
        irow = np.tile(np.arange(self.shape[0]), (self.shape[1], 1)).T
        irow -= self.scp_index[0]
        xrow = irow * self.grid_mult[0]

        icol = np.tile(np.arange(self.shape[1]), (self.shape[0], 1))
        icol -= self.scp_index[1]
        ycol = icol * self.grid_mult[1]
        return xrow, ycol

    def load_data(self):
        """Load data from the UMBRA SICD file."""
        reader = SICDReader(str(self.file_path))
        data = reader[:, :]
        return data

    def load_corrected_data(self, correction='sigma0'):
        """Load data with a specified correction"""
        if correction == 'sigma0':
            coeff = self.sigma0_coeff
        elif correction == 'beta0':
            coeff = self.beta0_coeff
        else:
            raise ValueError(f'Unknown correction type: {correction}')
        xrow, ycol = self.get_xrow_ycol()
        scale_factor = polyval2d(xrow, ycol, coeff)
        data = self.load_data()
        power = data.real**2 + data.imag**2
        corrected = power * scale_factor
        return corrected


def prep_umbra(granule_path: Path, work_dir: Optional[Path] = None) -> Path:
    """Prepare data for burst-based processing.

    Args:
        granule_path: Path to the UMBRA SICD file
        work_dir: Working directory for processing
    """
    if work_dir is None:
        work_dir = Path.cwd()
    reader = SICDReader(str(granule_path))
    sicd = reader.get_sicds_as_tuple()[0]
    umbra_sicd = UmbraSICD.from_sarpy_sicd(sicd, granule_path)

    dem_path = work_dir / 'dem.tif'
    dem.download_opera_dem_for_footprint(dem_path, umbra_sicd.footprint)
    return umbra_sicd, dem_path
