import numpy as np
from numpy.polynomial.polynomial import polyval2d
from osgeo import gdal

from multirtc.prep_umbra import UmbraSICD


def get_xrow_ycol(umbra_sicd: UmbraSICD) -> np.ndarray:
    """Calculate xrow and ycol for the umbra_sicd."""
    shadows_down_n_rows = umbra_sicd.shape[0]
    shadows_down_n_cols = umbra_sicd.shape[1]

    irow = np.tile(np.arange(shadows_down_n_cols), (shadows_down_n_rows, 1)).T
    irow -= umbra_sicd.scp_index[0]
    xrow = irow * umbra_sicd.azimuth_step

    icol = np.tile(np.arange(shadows_down_n_rows), (shadows_down_n_cols, 1))
    icol -= umbra_sicd.scp_index[1]
    ycol = icol * umbra_sicd.range_step
    return xrow, ycol


def save_as_beta0(umbra_sicd: UmbraSICD, output_dir) -> np.ndarray:
    """Save Umbra SICD as a beta0 geotiff file."""
    output_path = output_dir / f'{umbra_sicd.id}_beta0.tif'
    if output_path.exists():
        return output_path

    # FIXME: ISCE3 requires the complex-valued beta0. Not sure if the SICD sqrt(beta0 scaling factor) is valid for
    # complex-valued data.
    data = umbra_sicd.load_data()
    xrow, ycol = get_xrow_ycol(umbra_sicd)
    beta0_scale_factor = polyval2d(xrow, ycol, umbra_sicd.beta0_coeff)
    data = (data * np.sqrt(beta0_scale_factor)).T
    driver = gdal.GetDriverByName('GTiff')
    length, width = data.shape
    out_ds = driver.Create(str(output_path), width, length, 1, gdal.GDT_CFloat32)
    out_ds.GetRasterBand(1).WriteArray(data)
    out_ds.FlushCache()
    return output_path
