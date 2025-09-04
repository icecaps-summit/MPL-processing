# %%
import numpy as np
import pandas as pd

import sys
import subprocess
import gzip
from glob import glob
import xarray as xr

import mpl

from socket import gethostname
hostname = gethostname()
if 'jasmin' in hostname:
    din  = '/gws/ssde/j25b/icecaps/ICECAPSarchive/mpl/'
    dout = '/home/users/vonw/data/mpl/'
elif ('nuia' in hostname) or ('ncas' in hostname):
    din  = '/Users/vonw/data/icecaps/mpl/'
    dout = '/Users/vonw/data/icecaps/mpl/'
else:
    print('Unknown host!')
    exit()

# %%
def mplgz2ARM(date):
    
    yyyy = date.strftime('%Y')
    yyyymmdd = date.strftime('%Y%m%d')

    # ....Store filenames for one day
    fns = sorted(glob(din + 'raw/' + yyyy + '/' + yyyymmdd + '*.mpl.gz'))
    
    # ....Decode each hour
    data = []
    for fn in fns:
        #print(fn)
        # Define the path to your .gz file and the desired output file
        try:
            # Open the .gz file in binary read mode ('rb')
            with gzip.open(fn, 'rb') as f_in:
                # Open the output file in binary write mode ('wb')
                with open(fn[:-3], 'wb') as f_out:
                    # Read chunks from the compressed file and write them to the output file
                    # This is efficient for large files as it avoids loading the entire content into memory
                    for chunk in f_in:
                        f_out.write(chunk)
            #print(f"Successfully decompressed '{fn}' to '{fn[:-3]}'.")
            mpl108 = mpl.MPL(fn[:-3])
            mpl108.readData()
            data.append(mpl108.to_xarray())
            
        except FileNotFoundError:
            print(f"Error: The file '{fn}' was not found.")
        except Exception as e:
            print(f"An error occurred during decompression: {e}")
    
    # ....Combine hours into a dataset
    ds = xr.concat(data, dim='time')
    
    # ....Remove unnecessary variables
    ds = ds.drop_vars(['unitNumber',
                       'version',
                        'year',
                        'month',
                        'day',
                        'hours',
                        'minutes',
                        'seconds',
                        'Time',
                        'numberChannels',
                        'numberBins',
                        'binTime',
                        'rangeCalibration',
                        'numberDataBins',
                        'scanScenarioFlag',
                        'numberBackgrdBins',
                        'azimuthAngle',
                        'elevationAngle',
                        'compassDegrees',
                        'polarizationV0',
                        'polarizationV1',
                        'aToDdataBadFlag',
                        'dataFileVersion',
                        'mcsMode',
                        'firstDataBin',
                        'systemType',
                        'syncPulsePerSec',
                        'firstBackgrdBin',
                        'secondaryHdrSize',
                        'gpsLatitude',
                        'gpsLongitude',
                    ])

    # ....Rename existing variables to ARM variable names
    ds = ds.rename({'shotsSum': 'nshots',
               'triggerFrequency': 'rep_rate',
               'energyMonitor': 'energy',
               'detectorTemp': 'temp_detector',
               'telescopeTemp': 'temp_telescope',
               'laserTemp': 'temp_laser',
               'backgroundAverage': 'mn_background_1',
               'backgroundStdDev': 'std_background_1',
               'backgrdAverage2': 'mn_background_2',
               'backgrdStdDev2': 'std_background_2',
               'cloudBaseHeight': 'initial_cbh',
               })

    # ....Create ARM variables
    ds['base_time']   = int(ds.time[0])/1e9
    ds['base_time']   = ds['base_time'].astype('int32')
    time_offset       = (ds.time - ds.time[0]).values / 1e9
    ds['time_offset'] = xr.DataArray(time_offset, dims=['time']).astype(float)
    hours             = [pd.Timestamp(int(ds.base_time.values), unit='s') - date + pd.Timedelta(time, unit='s') for time in ds.time_offset.values]
    hour              = np.array([hour.seconds/3600 for hour in hours])
    ds['hour']        = xr.DataArray(hour, dims=['time']).astype('float32')
    ds['lat']         = 72.59622
    ds['lon']         = -38.42197
    ds['alt']         = 3200.          # Used this value to be consistent throughout entire ICECAPS period
    ds['nshots']      = ds['nshots'].astype('int32')
    for v in ['height', 'energy', 'temp_detector', 'temp_telescope', 'temp_laser', 'backscatter_1', 'backscatter_2', 'lat', 'lon', 'alt']:
        ds[v] = ds[v].astype('float32')

    # ....Add ARM global variables
    git = subprocess.run(['git', 'rev-parse', 'main'], capture_output=True)
    ds.attrs = {
        'Date_created': pd.Timestamp.utcnow().strftime('%a %b %d %H:%M:%S %Y UTC'),
        'Ingest_version': 'Git commit hash (SHA): ' + git.stdout.decode('utf-8'.rstrip()),
        'GitHub_repository': ''
        'comment': 'DOE Atmospheric Radiation Measurement (ARM) Micropulse Lidar (MPL) deployed to Summit, Greenland, as part of the NSF-funded ICECAPS project',
        'Author': 'Von P. Walden, Washington State University, v.walden@wsu.edu',
        'instrument_serial_number': '108',
        'instrument_version': '413',
        'backscatter_comment': 'See Flynn et al. 2007 Optics Express paper for details on how to interpret the two backscatter profiles'
    }
    # ....Add ARM varilable attributes
    ds['base_time'].attrs = {
        'long_name': 'Base time in Epoch',
        'units':     'seconds since 1970-1-1 0:00:00 0:00',
        }
    ds['time_offset'].attrs = {
        'long_name': 'Time offset from base_time',
        'units':     'seconds',
        }
    ds['hour'].attrs = {
        'long_name': 'Hour of the day',
        'units': 'UTC'
    }
    ds['height'].attrs = {
        'long_name': 'height',
        'units': 'km AGL'
    }
    ds['nshots'].attrs = {
        'long_name': 'number of laser shots',
        'units': 'unitless'
    }
    ds['rep_rate'].attrs = {
        'long_name': 'laser pulse repetition frequency',
        'units': 'Hz'
    }
    ds['energy'].attrs = {
        'long_name': 'laser energy',
        'units': 'microJoules'
    }
    ds['temp_detector'].attrs = {
        'long_name': 'detector temperature',
        'units': 'C'
    }
    ds['temp_telescope'].attrs = {
        'long_name': 'telescope temperature',
        'units': 'C'
    }
    ds['temp_laser'].attrs = {
        'long_name': 'laser temperature',
        'units': 'C'
    }
    ds['mn_background_1'].attrs = {
        'long_name': 'ean background in channel 1',
        'units': 'counts / microsecond'
    }
    ds['std_background_1'].attrs = {
        'long_name': 'standard deviation of the background in channel 1',
        'units': 'counts / microsecond'
    }
    ds['mn_background_2'].attrs = {
        'long_name': 'mean background in channel 2',
        'units': 'counts / microsecond'
    }
    ds['std_background_2'].attrs = {
        'long_name': 'standard deviation of the background in channel 2',
        'units': 'counts / microsecond'
    }
    ds['initial_cbh'].attrs = {
        'long_name': 'initial cloud base height from MPL software',
        'units': 'km AGL'
    }
    ds['backscatter_1'].attrs = {
        'long_name': 'attenuated backscatter in channel 1',
        'units': 'counts / microsecond',
        'channel_interpretation': 'This is the linear cross-polarization channel.  It is sensitive to the depolarized backscatter from the atmosphere',
        'comment': 'This field literally contains the counts detected by the detector for each range bin.  No corrections of any kind have been applied to this field.  In order to make proper use of the data, one should correct for detector non-linearity, subtract the afterpulse, subtract background counts, apply a range-squared correction, and correct for optical overlap and collimation effects'
    }
    ds['backscatter_2'].attrs = {
        'long_name': 'attenuated backscatter in channel 2',
        'units': 'counts / microsecond',
        'channel_interpretation': 'This is the circular polarization channel.  It is sensitive to the unpolarized backscatter from the atmosphere',
        'comment': 'This field literally contains the counts detected by the detector for each range bin.  No corrections of any kind have been applied to this field.  In order to make proper use of the data, one should correct for detector non-linearity, subtract the afterpulse, subtract background counts, apply a range-squared correction, and correct for optical overlap and collimation effects'
    }
    ds['lat'].attrs = {
        'long_name': 'north latitude',
        'units': 'deg'
    }
    ds['lon'].attrs = {
        'long_name': 'east longitude',
        'units': 'deg'
    }
    ds['alt'].attrs = {
        'long_name': 'altitude',
        'units': 'm MSL'
    }
    
    ds.to_netcdf(dout + 'ingested/smtmplpolX1.a1.' + yyyymmdd + '.000000.cdf')

    return

# %%
dates = pd.date_range('2017-01-21', '2017-01-22')
#dates = pd.date_range('2010-06-01', '2010-06-01')

for date in dates:
    try:
        print('Processing: ', date.strftime('%Y%m%d'))
        mplgz2ARM(date)
    except RuntimeWarning:
        print('SKIPPING:', date.strftime('%Y%m%d'))

# %%
