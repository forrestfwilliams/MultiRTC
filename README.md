# MultiRTC

A python library for creating ISCE3-based RTCs for multiple SAR data sources

> [!WARNING]
> This package is still in early development. Users are encouraged to not use this package in production or other critical contexts until the v1.0.0 release.

> [!IMPORTANT]
> All credit for this library's RTC algorithm goes to Gustavo Shiroma and the JPL [OPERA](https://www.jpl.nasa.gov/go/opera/about-opera/) and [ISCE3](https://github.com/isce-framework/isce3) teams. This package merely allows others to use their algorithm with a wider set of SAR data sources. The RTC algorithm utilized by this package is described in [Shiroma et al., 2023](https://doi.org/10.1109/TGRS.2022.3147472).

## Dataset Support
MultiRTC allows users to create RTC products from SLC data for multiple SAR sensor platforms, and provides utilities for assessing the resulting products. All utilities can be accessed via CLI pattern `multirtc SUBCOMMAND ARGS`, with the primary subcommand `multirtc rtc`.

Below is a list of relevant SAR data sources and their support status:

| Mission    | File Format | Image Mode        | Image Grid Type | Status      |
|------------|-------------|-------------------------------------|-------------|
| Sentinel-1 | SAFE        | Burst IW          | Range Doppler   | Supported   |
| Sentinel-1 | SAFE        | Full-frame IW     | Range Doppler   | Unsupported |
| Sentinel-1 | SAFE        | Burst EW          | Range Doppler   | Unsupported |
| Sentinel-1 | SAFE        | Full-frame EW     | Range Doppler   | Unsupported |
| Capella    | SICD        | Spotlight         | Polar           | Unsupported |
| Capella    | SICD        | Sliding Spotlight | Range Doppler   | Supported   |
| Capella    | SICD        | Stripmap          | Range Doppler   | Supported   | 
| Iceye      | SICD        | Dwell             | Range Doppler   | Supported   | 
| Iceye      | SICD        | Spotlight         | Range Doppler   | Supported   | 
| Iceye      | SICD        | Sliding Spotlight | Range Doppler   | Supported   |
| Iceye      | SICD        | Stripmap          | Range Doppler   | Supported   |
| Iceye      | SICD        | Scan              | Range Doppler   | Supported   |
| Umbra      | SICD        | Dwell             | Polar           | Supported^* |
| Umbra      | SICD        | Spotlight         | Polar           | Supported^* |

I have done my best to accurately reflect the support status of each SAR image type, but please let me know if I have made any mistakes. Note that some commercial datasets used to use polar instead of range doppler image grids for specific images modes. This table is based on the image grid types currently being used.

^*Polar image grid support is implemented via the [approach detailed by Piyush Agram](https://arxiv.org/abs/2503.07889v1) in his recent technical note. I have implemented his method in a fork of the main ISCE3 repo, which you can view [here](https://github.com/forrestfwilliams/isce3/tree/pfa). The long-term plan is to merge this into the main ISCE3 repo but until that is complete, polar grid support is only available via this project's `pfa`-suffixed docker containers. See the running via docker section for more details.

## Usage

To create an RTC, use the `multirtc` CLI entrypoint using the following pattern:

```bash
multirtc rtc PLATFORM SLC-GRANULE --resolution RESOLUTION --work-dir WORK-DIR
```
Where `PLATFORM` is the name of the satellite platform (currently `S1`, `CAPELLA`, `ICEYE` or `UMBRA`), `SLC-GRANULE` is the name of the SLC granule, `RESOLUTION` is the desired output resolution of the RTC image in meters, and `WORK-DIR` is the name of the working directory to perform processing in. Inputs such as the SLC data, DEM, and external orbit information are stored in `WORK-DIR/input`, while the RTC image and associated outputs are stored in `WORK-DIR/output` once processing is complete. SLC data that is available in the [Alaska Satellite Facility's data archive](https://search.asf.alaska.edu/#/?maxResults=250) (such as Sentinel-1 Burst SLCs) will be automatically downloaded to the input directory, but data not available in this archive (commercial datasets) are required to be staged in the input directory prior to processing.

Output RTC pixel values represent gamma0 power.

### Running via Docker
In addition to the main python interface, I've also provided an experimental docker container that contains full support for polar grid format SICD data. Encapsulating this functionality in a docker container is ncessary for now because it requires re-compiling a development version of ISCE3. The docker container can be run using a similar interface, with exception of needing to pass your EarthData credentials and the need to pass a mounted volume with an `input` and `output` directory inside:

```bash
docker run -it --rm \
    -e EARTHDATA_USERNAME=[YOUR_USERNAME_HERE] \
    -e EARTHDATA_PASSWORD=[YOUR_PASSWORD_HERE] \
    -v ~/local_path/project1:/home/conda/project1 \
    ghcr.io/forrestfwilliams/multirtc:0.4.0.pfa \
    rtc PLATFORM SLC-GRANULE --resolution RESOLUTION
```
The local `project1` directory can be a name of your choosing and should have the structure:
```
project1/
    |--input/
        |--input.slc (if needed)
    |--output/
```

### DEM options
Currently, only the OPERA DEM is supported. This is a global Height Above Ellipsoid DEM sourced from the [COP-30 DEM](https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.3). In the future, we hope to support a wider variety of automatically retrieved and user provided DEMs.

## Calibration & Validation Subcommands
MultiRTC includes three calibration and validation (cal/val) subcommands for assessing the geometric and radiometric quality of SAR products. These tools are useful for analyzing geolocation, co-registration, and impulse response performance.

### `ale` Absolute Location Error
Quantifies the geolocation accuracy of a SAR image by comparing known corner reflectors at the Rosamond, California site with their positions in the geocoded image.

Usage:
```bash
multirtc ale FILEPATH DATE AZMANGLE PROJECT --basedir BASEDIR
```
See `multirtc ale --help` for descriptions of each argument.

### `rle` Relative Location Error
Measures the relative alignment of overlapping geocoded SAR images by measuring the offsets between each 1024x1024 pixel chunk of the images.

Usage:
```bash
multirtc rle REFPATH SECPATH PROJECT --basedir BASEDIR
```
See `multirtc rle --help` for descriptions of each argument.

### `pt` Point Target Analysis
Evaluates the impulse response of corner reflector at the Rosamond, California site in the SAR image, including resolution, peak to side-lobe ratio (PSLR), and integrated side-lobe ratio (ISLR).

Usage:
```bash
multirtc pt FILEPATH DATE AZMANGLE PROJECT --basedir BASEDIR
```
See `multirtc pt --help` for descriptions of each argument.

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
