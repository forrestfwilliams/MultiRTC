import numpy as np
from numpy.polynomial.polynomial import polyval2d
from osgeo import gdal

from multirtc.prep_umbra import UmbraSICD


def get_xrow_ycol(umbra_sicd: UmbraSICD) -> np.ndarray:
    """Calculate xrow and ycol for the umbra_sicd."""
    irow = np.tile(np.arange(umbra_sicd.shape[0]), (umbra_sicd.shape[1], 1)).T
    irow -= umbra_sicd.scp_index[0]
    xrow = irow * umbra_sicd.azimuth_step

    icol = np.tile(np.arange(umbra_sicd.shape[1]), (umbra_sicd.shape[0], 1))
    icol -= umbra_sicd.scp_index[1]
    ycol = icol * umbra_sicd.range_step
    return xrow, ycol


def save_as_beta0(umbra_sicd: UmbraSICD, output_dir) -> np.ndarray:
    """Save Umbra SICD as a beta0 geotiff file."""
    output_path = output_dir / f'{umbra_sicd.id}_beta0.tif'
    if output_path.exists():
        return output_path

    data_base = umbra_sicd.load_data()
    data_power = data_base.real**2 + data_base.imag**2

    xrow, ycol = get_xrow_ycol(umbra_sicd)
    beta0_scale_factor = polyval2d(xrow, ycol, umbra_sicd.beta0_coeff)

    data_beta0 = data_power * beta0_scale_factor

    driver = gdal.GetDriverByName('GTiff')
    length, width = data_power.shape
    out_ds = driver.Create(str(output_path), width, length, 1, gdal.GDT_Float32)
    out_ds.GetRasterBand(1).WriteArray(data_beta0)
    out_ds.FlushCache()
    return output_path
