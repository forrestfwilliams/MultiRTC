from pathlib import Path
from typing import Union

import requests


def download_file(
    url: str,
    download_path: Union[Path, str] = '.',
    chunk_size=10 * (2**20),
) -> Path:
    """Download a file without authentication.

    Args:
        url: URL of the file to download
        download_path: Path to save the downloaded file to
        chunk_size: Size to chunk the download into

    Returns:
        download_path: The path to the downloaded file
    """
    session = requests.Session()

    with session.get(url, stream=True) as s:
        s.raise_for_status()
        with open(download_path, 'wb') as f:
            for chunk in s.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
    session.close()
