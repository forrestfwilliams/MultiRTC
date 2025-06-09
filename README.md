# MultiRTC

A python library for creating ISCE3-based RTCs for multiple SAR data sources

> [!WARNING]
> This package is still in early development. Users are encouraged to not use this package in production or other critical contexts until the v1.0.0 release.

> [!IMPORTANT]
> All credit for this library's RTC algorithm goes to Gustavo Shiroma and the JPL [OPERA](https://www.jpl.nasa.gov/go/opera/about-opera/) and [ISCE3](https://github.com/isce-framework/isce3) teams. This package merely allows others to use their algorithm with a wider set of SAR data sources. The RTC algorithm utilized by this package is described in [Shiroma et al., 2023](https://doi.org/10.1109/TGRS.2022.3147472).

## Usage
MultiRTC allows users to create RTC products from SLC data for multiple SAR sensor platforms. Currently this list includes:

Full RTC:
- [Sentinel-1 Burst SLCs](https://www.earthdata.nasa.gov/data/catalog/alaska-satellite-facility-distributed-active-archive-center-sentinel-1-bursts-version)
- [Capella SICD SLCs](https://www.capellaspace.com/earth-observation/data)
- [ICEYE SICD SLCs](https://sar.iceye.com/5.0/productFormats/slc/)

Geocode Only:
- [UMBRA SICD SLCs](https://help.umbra.space/product-guide/umbra-products/umbra-product-specifications)

To create an RTC, use the `multirtc` CLI entrypoint using the following pattern:

```bash
multirtc PLATFORM SLC-GRANULE --resolution RESOLUTION --work-dir WORK-DIR
```
Where `PLATFORM` is the name of the satellite platform (currently `S1`, `CAPELLA`, `ICEYE` or `UMBRA`), `SLC-GRANULE` is the name of the SLC granule, `RESOLUTION` is the desired output resolution of the RTC image in meters, and `WORK-DIR` is the name of the working directory to perform processing in. Inputs such as the SLC data, DEM, and external orbit information are stored in `WORK-DIR/input`, while the RTC image and associated outputs are stored in `WORK-DIR/output` once processing is complete. SLC data that is available in the [Alaska Satellite Facility's data archive](https://search.asf.alaska.edu/#/?maxResults=250) (such as Sentinel-1 Burst SLCs) will be automatically downloaded to the input directory, but data not available in this archive (commercial datasets) are required to be staged in the input directory prior to processing.

Output RTC pixel values represent gamma0 power.

### Current Umbra Implementation
Currently, the Umbra processor only supports basic geocoding and not full RTC processing. ISCE3's RTC algorithm is only designed to work with Range Migration Algorithm (RMA) focused SLC products, but Umbra creates their data using the Polar Format Algorithm (PFA). Using an [approach detailed by Piyush Agram](https://arxiv.org/abs/2503.07889v1) to adapt RMA approaches to the PFA image geometry, we have developed a workflow to geocode an Umbra SLC but there is more work to be done to implement full RTC processing. Since full RTC is not yet implemented, Umbra geocoded pixel values represent sigma0 power.

### DEM options
Currently, only the OPERA DEM is supported. This is a global Height Above Ellipsoid DEM sourced from the [COP-30 DEM](https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.3). In the future, we hope to support a wider variety of automatically retrieved and user provided DEMs.

## When will support for [insert SAR provider here] products be added?
We're currently working on this package on a "best effort" basis with no specific timeline for any particular dataset. We would love to add support for every SAR dataset ASAP, but we only have so much time to devote to this package. If you want a particular dataset to be prioritized there are several things you can do:

- [Open an issue](https://github.com/forrestfwilliams/multirtc/issues/new) requesting support for your dataset and encourage others to like or comment on it.
- Provides links to example datasets over the Rosamond, California corner reflector site (Lat/Lon 34.799,-118.095) for performing cal/val.
- Reach out to us about funding the development required to add your dataset.

## Developer Setup
1. Ensure that conda is installed on your system (we recommend using [mambaforge](https://github.com/conda-forge/miniforge#mambaforge) to reduce setup times).
2. Download a local version of the `multirtc` repository (`git clone https://github.com/forrestfwilliams/multirtc.git`)
3. In the base directory for this project call `mamba env create -f environment.yml` to create your Python environment, then activate it (`mamba activate multirtc`)
4. Finally, install a development version of the package (`python -m pip install -e .`)

To run all commands in sequence use:
```bash
git clone https://github.com/forrestfwilliams/multirtc.git
cd multirtc
mamba env create -f environment.yml
mamba activate multirtc
python -m pip install -e .
```

## License
MultiRTC is licensed under the BSD-3-Clause license. See the LICENSE file for more details.

## Code of conduct
We strive to create a welcoming and inclusive community for all contributors to this project. As such, all contributors to this project are expected to adhere to our code of conduct.

Please see `CODE_OF_CONDUCT.md` for the full code of conduct text.

## Contributing
Contributions to this project plugin are welcome! If you would like to contribute, please submit a pull request on the GitHub repository.

## Contact Us
Want to talk about this project? We would love to hear from you!

Found a bug? Want to request a feature?
[open an issue](https://github.com/forrestfwilliams/multirtc/issues/new)
