from abc import ABC
from datetime import datetime
from pathlib import Path

import isce3
import numpy as np
from shapely.geometry import Point, Polygon


def to_isce_datetime(dt):
    if isinstance(dt, datetime):
        return isce3.core.DateTime(dt)
    elif isinstance(dt, np.datetime64):
        return isce3.core.DateTime(dt.item())
    else:
        raise ValueError(f'Unsupported datetime type: {type(dt)}. Expected datetime or np.datetime64.')


def from_isce_datetime(dt):
    return datetime.fromisoformat(dt.isoformat())


class SlcTemplate(ABC):
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
