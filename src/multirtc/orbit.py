import os
from datetime import datetime, timedelta
from pathlib import Path
from typing import Union

import boto3
from botocore import UNSIGNED
from botocore.client import Config
from hyp3lib.fetch import download_file


s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

# orbits bucket
BUCKET = 's1-orbits'

# date format used in file names
FMT = '%Y%m%dT%H%M%S'

# Orbital period of Sentinel-1 in seconds:
# 12 days * 86400.0 seconds/day, divided into 175 orbits
T_ORBIT_INT = (12 * 86400.0) / 175.0
PADDING_SHORT_INT = 60

# Temporal margin to apply to the start time of a frame
#  to make sure that the ascending node crossing is
#    included when choosing the orbit file
PADDING_SHORT = timedelta(seconds=PADDING_SHORT_INT)
MARGIN_START_TIME = timedelta(seconds=T_ORBIT_INT + PADDING_SHORT_INT)


def build_url(bucket: str, key: str) -> str:
    return f'https://{bucket}.s3.amazonaws.com/{key}'


def list_bucket(bucket: str, prefix: str) -> list[str]:
    paginator = s3.get_paginator('list_objects_v2')
    page_iterator = paginator.paginate(
        Bucket=bucket,
        Prefix=prefix,
    )
    keys = [item['Key'] for page in page_iterator for item in page.get('Contents', [])]
    keys.sort(reverse=True)
    return keys


def get_orbit_for_granule(granule: str, bucket: str, orbit_type: str) -> Union[str, None]:
    platform = granule[0:3]
    granule_start_date = datetime.strftime(datetime.strptime(granule[17:32], FMT) - MARGIN_START_TIME, FMT)
    granule_end_date = datetime.strftime(datetime.strptime(granule[33:48], FMT) + PADDING_SHORT, FMT)

    keys = list_bucket(bucket=bucket, prefix=f'{orbit_type}/{platform}')
    for key in keys:
        filename = os.path.basename(key)
        orbit_start_date = filename[42:57]
        orbit_end_date = filename[58:73]
        if orbit_start_date <= granule_start_date <= granule_end_date <= orbit_end_date:
            return key
    return None


def get_url(granule: str, bucket: str) -> Union[str, None]:
    for orbit_type in ['AUX_POEORB', 'AUX_RESORB']:
        key = get_orbit_for_granule(granule, bucket, orbit_type)
        if key:
            return build_url(bucket, key)
    return None


def get_orbit(scene: str, save_dir: Path):
    url = get_url(scene, BUCKET)
    orbit_path = save_dir / url.split('/')[-1]

    if orbit_path.exists():
        return orbit_path

    download_file(url, orbit_path.parent)
    return orbit_path
