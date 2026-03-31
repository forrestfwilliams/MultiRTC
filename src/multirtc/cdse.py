"""Download Sentinel-1 SLC products from the Copernicus Data Space Ecosystem (CDSE).

This module provides an alternative to ASF-based downloads via burst2safe.
It searches the CDSE OData catalog by the parent SLC scene name (derived from
the burst granule via ASF search) and downloads the full SLC zip.

CDSE credentials (username/password) can be provided via:
  - Environment variables: CDSE_USERNAME and CDSE_PASSWORD
  - The ~/.netrc file with machine: dataspace.copernicus.eu

References:
  - https://documentation.dataspace.copernicus.eu/APIs/OData.html
  - https://documentation.dataspace.copernicus.eu/APIs/Token.html
"""

import logging
import netrc
import os
import time
import zipfile
from pathlib import Path

import asf_search
import requests


logger = logging.getLogger(__name__)

CDSE_HOST = 'dataspace.copernicus.eu'
CDSE_TOKEN_URL = 'https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token'
CDSE_ODATA_URL = 'https://catalogue.dataspace.copernicus.eu/odata/v1/Products'
CDSE_DOWNLOAD_URL = 'https://download.dataspace.copernicus.eu/odata/v1/Products'

# Retry settings
MAX_RETRIES = 3
RETRY_START_WAIT = 10  # seconds
RETRY_INCREMENT = 10  # seconds


def get_cdse_credentials(
    username: str | None = None,
    password: str | None = None,
) -> tuple[str, str]:
    """Resolve CDSE credentials from arguments, environment, or ~/.netrc.

    Args:
        username: CDSE username. Falls back to CDSE_USERNAME env var, then ~/.netrc.
        password: CDSE password. Falls back to CDSE_PASSWORD env var, then ~/.netrc.

    Returns:
        Tuple of (username, password).
    """
    if username and password:
        return username, password

    env_user = os.getenv('CDSE_USERNAME')
    env_pass = os.getenv('CDSE_PASSWORD')
    if env_user and env_pass:
        return env_user, env_pass

    try:
        nrc = netrc.netrc()
        auth = nrc.authenticators(CDSE_HOST)
        if auth:
            return auth[0], auth[2]
    except (FileNotFoundError, netrc.NetrcParseError):
        pass

    raise ValueError(
        'CDSE credentials not found. Provide them via:\n'
        '  1. CDSE_USERNAME and CDSE_PASSWORD environment variables\n'
        '  2. ~/.netrc entry for machine dataspace.copernicus.eu\n'
        'Register for a free account at https://dataspace.copernicus.eu/'
    )


def ensure_cdse_credentials(
    username: str | None = None,
    password: str | None = None,
) -> None:
    """Ensure CDSE credentials are available in ~/.netrc.

    If credentials are provided via env vars but ~/.netrc does not
    contain an entry for CDSE, the entry will be appended to ~/.netrc.
    """
    if username is None:
        username = os.getenv('CDSE_USERNAME')
    if password is None:
        password = os.getenv('CDSE_PASSWORD')

    netrc_file = Path.home() / '.netrc'

    cdse_in_netrc = False
    if netrc_file.exists():
        try:
            nrc = netrc.netrc(netrc_file)
            if nrc.authenticators(CDSE_HOST):
                cdse_in_netrc = True
        except netrc.NetrcParseError:
            pass

    if username and password and not cdse_in_netrc:
        with open(netrc_file, 'a') as f:
            f.write(f'\nmachine {CDSE_HOST} login {username} password {password}\n')
        netrc_file.chmod(0o600)
    elif username and password and cdse_in_netrc:
        logging.info(f'CDSE credentials already present in {netrc_file}, skipping update.')

    get_cdse_credentials(username, password)


def get_cdse_access_token(username: str, password: str) -> str:
    """Obtain an access token from the CDSE identity provider.

    Args:
        username: CDSE username.
        password: CDSE password.

    Returns:
        Access token string.
    """
    data = {
        'grant_type': 'password',
        'username': username,
        'password': password,
        'client_id': 'cdse-public',
    }
    response = requests.post(CDSE_TOKEN_URL, data=data, timeout=60)
    response.raise_for_status()
    return response.json()['access_token']


def search_cdse_by_scene_name(scene_name: str) -> dict:
    """Search the CDSE OData catalog for a Sentinel-1 SLC by scene name.

    Args:
        scene_name: Sentinel-1 scene name (without .SAFE or .zip extension).

    Returns:
        Product entry from the CDSE OData response containing 'Id' and 'Name'.

    Raises:
        LookupError: If the product is not found on CDSE.
    """
    scene_name = scene_name.replace('.zip', '').replace('.SAFE', '')
    safe_name = f'{scene_name}.SAFE'
    query = f"{CDSE_ODATA_URL}?$filter=Name eq '{safe_name}'"

    response = requests.get(query, timeout=120)
    response.raise_for_status()
    results = response.json().get('value', [])

    if not results:
        raise LookupError(f"Product '{safe_name}' not found in CDSE catalog.")
    return results[0]


def burst_to_parent_slc(burst_granule: str) -> str:
    """Use ASF search to find the parent SLC scene name for a burst granule.

    Args:
        burst_granule: Burst granule name (e.g. S1_136231_IW2_20200604T022312_VV_7C85-BURST).

    Returns:
        Parent SLC scene name (e.g. S1A_IW_SLC__1SDV_20200604T022251_20200604T022318_032861_03CE65_7C85).
    """
    results = asf_search.granule_search([burst_granule])
    if not results:
        raise LookupError(f'Burst granule {burst_granule} not found in ASF archive.')

    url = results[0].properties['url']
    # URL format: https://sentinel1-burst.asf.alaska.edu/{PARENT_SLC}/IW{N}/{POL}/{idx}.tiff
    parent_slc = url.split('/')[3]
    logger.info(f'Mapped burst {burst_granule} to parent SLC {parent_slc}')
    return parent_slc


def download_slc_from_cdse(
    scene_name: str,
    output_dir: Path | str,
    max_retries: int = MAX_RETRIES,
) -> Path:
    """Download a Sentinel-1 SLC product from CDSE and extract the SAFE directory.

    Args:
        scene_name: Sentinel-1 scene name.
        output_dir: Directory to save and extract the downloaded product.
        max_retries: Number of download attempts before raising an error.

    Returns:
        Path to the extracted .SAFE directory.
    """
    output_dir = Path(output_dir).resolve()
    scene_name = scene_name.replace('.zip', '').replace('.SAFE', '')
    safe_dir = output_dir / f'{scene_name}.SAFE'

    # Skip download if SAFE already exists
    if safe_dir.exists():
        logger.info(f'SAFE directory already exists: {safe_dir}')
        return safe_dir

    # Get credentials and token
    cdse_user, cdse_pass = get_cdse_credentials()
    access_token = get_cdse_access_token(cdse_user, cdse_pass)

    # Search CDSE catalog
    product = search_cdse_by_scene_name(scene_name)
    product_id = product['Id']

    download_url_zip = f'{CDSE_DOWNLOAD_URL}({product_id})/$zip'
    download_url_value = f'{CDSE_DOWNLOAD_URL}({product_id})/$value'
    headers = {'Authorization': f'Bearer {access_token}'}

    out_zip = output_dir / f'{scene_name}.zip'

    def _do_download(url: str) -> None:
        """Perform the actual download from a given URL."""
        response = requests.get(url, headers=headers, stream=True, timeout=600)
        response.raise_for_status()

        with open(out_zip, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192 * 16):
                if chunk:
                    f.write(chunk)

        if out_zip.stat().st_size == 0:
            out_zip.unlink(missing_ok=True)
            raise requests.RequestException('Downloaded file is empty')

        logger.info(f'Downloaded {out_zip.name} from CDSE ({out_zip.stat().st_size / 1e6:.1f} MB)')

    last_exc: Exception | None = None
    for attempt in range(1, max_retries + 1):
        logger.info(f'CDSE download attempt #{attempt} for {scene_name}')
        try:
            try:
                _do_download(download_url_zip)
                break
            except requests.HTTPError as e:
                out_zip.unlink(missing_ok=True)
                if e.response is not None and e.response.status_code == 404:
                    logger.info('Compressed format not available, falling back to uncompressed...')
                    try:
                        _do_download(download_url_value)
                        break
                    except requests.RequestException:
                        out_zip.unlink(missing_ok=True)
                        raise
                raise
            except requests.RequestException:
                out_zip.unlink(missing_ok=True)
                raise
        except Exception as exc:
            last_exc = exc
            wait_time = RETRY_START_WAIT + RETRY_INCREMENT * (attempt - 1)
            if attempt < max_retries:
                logger.warning(f'Attempt #{attempt} failed: {exc}. Waiting {wait_time}s before retry...')
                time.sleep(wait_time)
    else:
        raise RuntimeError(
            f'Failed to download {scene_name} from CDSE after {max_retries} attempts'
        ) from last_exc

    # Extract the zip to get the SAFE directory
    logger.info(f'Extracting {out_zip.name}...')
    with zipfile.ZipFile(out_zip, 'r') as zf:
        zf.extractall(output_dir)
    out_zip.unlink()

    if not safe_dir.exists():
        # Some CDSE zips may have a different top-level name; find the .SAFE directory
        safe_dirs = list(output_dir.glob('*.SAFE'))
        if safe_dirs:
            safe_dir = safe_dirs[0]
        else:
            raise FileNotFoundError(f'Could not find extracted SAFE directory in {output_dir}')

    logger.info(f'Extracted SAFE directory: {safe_dir}')
    return safe_dir
