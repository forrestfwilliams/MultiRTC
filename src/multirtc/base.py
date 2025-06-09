from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path

import isce3
import numpy as np
from shapely.geometry import Point, Polygon


def to_isce_datetime(dt: datetime | np.datetime64) -> isce3.core.DateTime:
    if isinstance(dt, datetime):
        return isce3.core.DateTime(dt)
    elif isinstance(dt, np.datetime64):
        return isce3.core.DateTime(dt.item())
    else:
        raise ValueError(f'Unsupported datetime type: {type(dt)}. Expected datetime or np.datetime64.')


def from_isce_datetime(dt: isce3.core.DateTime) -> datetime:
    return datetime.fromisoformat(dt.isoformat())


def print_wkt(slc):
    radar_grid = slc.radar_grid
    dem = isce3.geometry.DEMInterpolator(slc.scp_hae)
    doppler = slc.doppler_centroid_grid
    wkt = isce3.geometry.get_geo_perimeter_wkt(
        grid=radar_grid, orbit=slc.orbit, doppler=doppler, dem=dem, points_per_edge=3
    )
    print(wkt)


class Slc(ABC):
    """Template class for SLC objects that defines a common interface and enforces required attributes."""

    required_attributes = {
        'id': str,
        'filepath': Path,
        'footprint': Polygon,
        'center': Point,
        'lookside': str,  # 'right' or 'left'
        'wavelength': float,
        'polarization': str,
        'shape': tuple,
        'range_pixel_spacing': float,
        'reference_time': datetime,
        'sensing_start': float,
        'prf': float,
        'supports_rtc': bool,
        'supports_bistatic_delay': bool,
        'supports_static_tropo': bool,
        'orbit': object,  # Replace with actual orbit type
        'radar_grid': object,  # Replace with actual radar grid type
        'doppler_centroid_grid': object,  # Replace with actual doppler centroid grid type
    }

    # I prefer this setup to enforce properties over forcing subclasses to have a bunch of @property statements
    def __init_subclass__(cls):
        super().__init_subclass__()
        original_init = cls.__init__

        def wrapped_init(self, *args, **kwargs):
            original_init(self, *args, **kwargs)
            for attr, expected_type in cls.required_attributes.items():
                if not hasattr(self, attr):
                    raise NotImplementedError(f'{cls.__name__} must define self.{attr}')
                if not isinstance(getattr(self, attr), expected_type):
                    raise TypeError(
                        f'{cls.__name__}.{attr} must be of type {expected_type.__name__},'
                        f'got {type(getattr(self, attr)).__name__}'
                    )

        cls.__init__ = wrapped_init

    @abstractmethod
    def create_geogrid(self, spacing_meters: int) -> isce3.product.GeoGridParameters:
        """
        Create a geogrid for the SLC object with the specified resolution.

        Args:
            spacing_meters: Pixel spacing in meters for the geogrid.

        Returns:
            The geogrid parameters for the SLC object.
        """
        pass
