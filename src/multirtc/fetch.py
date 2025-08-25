"""Utilities for fetching things from external endpoints. Vendored from hyp3lib"""

import logging
from email.message import Message
from os.path import basename
from pathlib import Path
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


def _get_download_path(url: str, content_disposition: str | None = None, directory: Path | str = '.'):
    filename = None
    if content_disposition is not None:
        message = Message()
        message['content-type'] = content_disposition
        filename = message.get_param('filename')
    if not filename:
        filename = basename(urlparse(url).path)
    if not filename:
        raise ValueError(f'could not determine download path for: {url}')
    assert isinstance(filename, str)
    return Path(directory) / filename


def download_file(
    url: str,
    directory: Path | str = '.',
    chunk_size=None,
    retries=2,
    backoff_factor=1,
    auth: tuple[str, str] | None = None,
    token: str | None = None,
) -> str:
    """Download a file

    Args:
        url: URL of the file to download
        directory: Directory location to place files into
        chunk_size: Size to chunk the download into
        retries: Number of retries to attempt
        backoff_factor: Factor for calculating time between retries
        auth: Username and password for HTTP Basic Auth
        token: Token for HTTP Bearer authentication

    Returns:
        download_path: The path to the downloaded file
    """
    logging.info(f'Downloading {url}')

    session = requests.Session()
    session.auth = auth
    if token:
        session.headers.update({'Authorization': f'Bearer {token}'})

    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    session.mount('https://', HTTPAdapter(max_retries=retry_strategy))
    session.mount('http://', HTTPAdapter(max_retries=retry_strategy))

    with session.get(url, stream=True) as s:
        download_path = _get_download_path(s.url, s.headers.get('content-disposition'), directory)
        s.raise_for_status()
        with open(download_path, 'wb') as f:
            for chunk in s.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    session.close()

    return str(download_path)
