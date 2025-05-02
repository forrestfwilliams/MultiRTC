# MultiRTC

A python library for creating ISCE3-based RTCs from Sentinel-1 Burst and Umbra SICD SLC data.

**ALL CREDIT FOR THIS LIBRARY'S RTC ALGORITHM GOES TO GUSTAVO SHIROMA AND THE [JPL OPERA TEAM](https://www.jpl.nasa.gov/go/opera). THIS PLUGIN MERELY ALLOWS OTHERS TO USE THEIR ALGORITHM WITH MULTIPLE SENSORS.**

## Usage
MultiRTC allows users to create RTC products from SLC data for multiple SAR sensor platforms. Currently this list includes:

- [Sentinel-1 burst SLCs](https://www.earthdata.nasa.gov/data/catalog/alaska-satellite-facility-distributed-active-archive-center-sentinel-1-bursts-version)
- [UMBRA SICD SLCs](https://help.umbra.space/product-guide/umbra-products/umbra-product-specifications)

To create an RTC, use the `multirtc` CLI entrypoint using the following pattern:

```bash
multirtc PLATFORM SLC-GRANULE --resolution RESOLUTION --work-dir WORK-DIR
```
Where `PLATFORM` is the name of the satellite platform (currently `S1` or `UMBRA`), `SLC-GRANULE` is the name of the SLC granule, `RESOLUTION` is the desired output resolution of the RTC image in meters, and `WORK-DIR` is the name of the working directory to perform processing in. Inputs such as the SLC data, DEM, and external orbit information are stored in `WORK-DIR/input`, while the RTC image and associated outputs are stored in `WORK-DIR/output` once processing is complete. SLC data that is available in the [Alaska Satellite Facility's data archive](https://search.asf.alaska.edu/#/?maxResults=250) (such as Sentinel-1 burst SLCs) will be automatically downloaded to the input directory, but data not available in this archive (Umbra SICD SLCs) are required to be staged in the input directory prior to processing.

Output RTC products are in Gamma0 radiometry.

### Current Umbra Implementation
Currently, the Umbra processor only supports basic geocoding and not full RTC processing. ISCE3's RTC algorithm is only designed to work with Range Migration Algorithm (RMA) focused SLC products, but Umbra creates their data using the Polar Format Algorithm (PFA). Using an [approach detailed by Piyush Agram](https://arxiv.org/abs/2503.07889v1) to adapt RMA approaches to the PFA image geometry, we have developed a workflow to geocode an Umbra SLC but there is more work to be done to implement full RTC processing. Since full RTC is not yet implemented, Umbra geocoded products are in Sigma0 radiometry.

### DEM options
Currently, only the NISAR DEM is supported. This is a roughly global Height Above Ellipsoid DEM sourced from the [COP-30 DEM](https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.3). In the future, we hope to support a wider variety of automatically retrieved and user provided DEMs.

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
