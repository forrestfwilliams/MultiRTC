import datetime
import os
from pathlib import Path

import isce3
import numpy as np
from isce3.stripmap.readers.l1.ALOS2.CEOS import ImageFile, LeaderFile
from isce3.stripmap.readers.l1.ALOS2.CEOS.SignalDataRecordType import SignalDataRecordType
from shapely import wkt

# from multirtc import define_geogrid
from multirtc.base import Slc, to_isce_datetime


class ScanImageFile(ImageFile.ImageFile):
    def readNextLine(self):
        # Create record type with information from description
        bytesperpixel = self.description.NumberOfBytesPerDataGroup // self.description.NumberOfSamplesPerDataGroup
        pixels = self.description.NumberOfBytesOfSARDataPerRecord // bytesperpixel
        record = SignalDataRecordType(pixels=pixels, bytesperpixel=bytesperpixel)
        # Read from file
        record.fromfile(self.fid)
        self.counter = self.counter + 1
        # Sensor specific validators
        assert record.RecordSequenceNumber == (self.counter + 1)
        assert record.FirstRecordType == 50
        assert record.RecordTypeCode == 10
        assert record.SecondRecordSubType == 18
        assert record.ThirdRecordSubType == 20
        assert record.ActualCountOfLeftFillPixels == 0
        assert record.SARChannelCode == 0
        # assert record.ScanIDForScanSAR == 0
        assert record.OnboardRangeCompressedFlag == 0
        assert record.PulseTypeDesignator == 0
        return record


def get_alos_filenames(product_path: Path):
    filenames = {}
    # First look for the leader file
    flist = list(product_path.glob('LED-ALOS2*1.1__*'))
    if len(flist) == 0:
        raise ValueError(f'No leader files found in folder {product_path}')
    elif len(flist) > 1:
        raise ValueError(f'Multiple leader files in folder {product_path}')

    filenames['leaderfile'] = flist[0]
    pattern = os.path.basename(flist[0])[4:]

    # FIXME add acutal list
    for pol in ['HH']:
        flist = sorted(list(product_path.glob(f'IMG-{pol}-{pattern}*')))
        assert len(flist) > 0, f'No image files found for polarization {pol} in folder {product_path}'
        filenames[pol] = flist

    filenames['defaulth5'] = f'{pattern}.h5'
    return filenames


def parse_leader_file(filenames):
    """
    Parse leader file and check values against polarizations.
    """
    try:
        ldr = LeaderFile.LeaderFile(filenames['leaderfile'])
    except AssertionError as msg:
        print(msg)
        raise AssertionError(f'Error parsing ALOS-2 L1.1 leader file: {filenames["leaderfile"]}')
    return ldr


class AL2ScanSlc(Slc):
    def __init__(self, product_path: Path, polarization='HH', swath=1):
        files = get_alos_filenames(product_path)
        leader = parse_leader_file(files)
        image = ScanImageFile(files[polarization][swath - 1])
        firstrec = image.readNextLine()
        bytesperpixel = image.description.NumberOfBytesPerDataGroup // image.description.NumberOfSamplesPerDataGroup
        width = (
            image.description.NumberOfBytesOfSARDataPerRecord // bytesperpixel
        ) // image.description.NumberOfSamplesPerDataGroup
        length = image.description.NumberOfSARDataRecords
        sensing_start_dt = datetime.datetime(firstrec.SensorAcquisitionYear, 1, 1) + datetime.timedelta(
            days=int(firstrec.SensorAcquisitionDayOfYear - 1), seconds=firstrec.SensorAcquisitionusecsOfDay * 1e-6
        )
        reference_time = self.get_refernce_time(leader)
        self.id = product_path.name
        self.filepath = product_path
        assert leader.summary.SensorIDAndMode[7] in ['L', 'R']
        self.lookside = 'right' if leader.summary.SensorIDAndMode[7] == 'R' else 'left'
        self.wavelength = leader.summary.RadarWavelengthInm
        self.polarization = polarization
        self.shape = (length, width)
        self.starting_range = firstrec.SlantRangeToFirstSampleInm
        self.range_pixel_spacing = isce3.core.speed_of_light / (2 * leader.summary.SamplingRateInMHz * 1.0e6)
        self.reference_time = to_isce_datetime(self.get_refernce_time(leader))
        self.sensing_start = (sensing_start_dt - reference_time).total_seconds()
        self.prf = firstrec.PRFInmHz * 1.0e-3
        self.supports_rtc = True
        self.supports_bistatic_delay = False
        self.supports_static_tropo = False
        self.orbit = self.get_orbit(leader)
        self.radar_grid = self.get_radar_grid()
        # FIXME: not sure if the data is truly zero-doppler
        self.doppler_centroid_grid = isce3.core.LUT2d()
        # self.doppler_centroid_grid = self.get_doppler_centroid_grid(leader)
        self.footprint = self.get_footprint()
        self.center = self.footprint.centroid

    def get_refernce_time(self, leader: LeaderFile.LeaderFile) -> datetime.datetime:
        hdr = leader.platformPosition.header
        return datetime.datetime(hdr.YearOfDataPoint, hdr.MonthOfDataPoint, hdr.DayOfDataPoint)

    def get_orbit(self, leader: LeaderFile.LeaderFile) -> isce3.core.Orbit:
        hdr = leader.platformPosition.header
        t0, dt = hdr.SecondsOfDay, hdr.TimeIntervalBetweenDataPointsInSec
        times, svs = [], []
        for i, sv in enumerate(leader.platformPosition.statevectors):
            t = t0 + i * dt * 1.0
            times.append(t)
            timestamp = self.reference_time + isce3.core.TimeDelta(seconds=t)
            svs.append(
                isce3.core.StateVector(
                    datetime=timestamp,
                    position=[sv.PositionXInm, sv.PositionYInm, sv.PositionZInm],
                    velocity=[sv.VelocityXInmpers, sv.VelocityYInmpers, sv.VelocityZInmpers],
                )
            )
        return isce3.core.Orbit(svs, self.reference_time, type='DOE')

    def get_radar_grid(self) -> isce3.product.RadarGridParameters:
        """Define the radar grid parameters for the SLC.

        Returns:
            An instance of isce3.product.RadarGridParameters representing the radar grid.
        """
        radar_grid = isce3.product.RadarGridParameters(
            sensing_start=self.sensing_start,
            wavelength=self.wavelength,
            prf=self.prf,
            starting_range=self.starting_range,
            range_pixel_spacing=self.range_pixel_spacing,
            lookside=isce3.core.LookSide.Right if self.lookside == 'right' else isce3.core.LookSide.Left,
            length=self.shape[0],
            width=self.shape[1],
            ref_epoch=self.reference_time,
        )
        return radar_grid

    def get_doppler_centroid_grid(self, leader):
        doppler_coeffs = [
            leader.summary.DopplerCenterFrequencyConstantTerm,
            leader.summary.DopplerCenterFrequencyLinearTerm,
        ]
        rng = self.starting_range + np.arange(0, self.shape[1], 100) * self.range_pixel_spacing
        doppler = doppler_coeffs[0] + doppler_coeffs[1] * rng / 1000
        dfit = np.polyfit(np.arange(0, self.shape[1], 100), doppler, 1)
        doppler_coeffs_rbin = [dfit[1], dfit[0], 0.0, 0.0]
        return doppler_coeffs_rbin

    def get_footprint(self):
        radar_grid = self.radar_grid
        dem = isce3.geometry.DEMInterpolator(1000)
        doppler = self.doppler_centroid_grid
        wkt_str = isce3.geometry.get_geo_perimeter_wkt(
            grid=radar_grid, orbit=self.orbit, doppler=doppler, dem=dem, points_per_edge=3
        )
        footprint = wkt.loads(wkt_str)
        return footprint

    def create_geogrid(self, spacing_meters: int) -> isce3.product.GeoGridParameters:
        # return define_geogrid.generate_geogrids(self, spacing_meters, self.local_epsg)
        return None
