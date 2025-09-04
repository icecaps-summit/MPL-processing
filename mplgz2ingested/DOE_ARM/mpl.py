"""
Class definition for the ARM Micro Pulse Lidar (MPL).

@author: Von P. Walden, 13 October 2014
         (C code that was created by Rich Coulter and Dave Turner was helpful
          when this Python code was written.)
          
          Usage:
          		>>> import mpl
          		>>> mpl108 = mpl.MPL('/Users/vonw/data/N-ICE/mpl/201411230000.mpl')
          		>>> mpl108.readData()
          		>>> dir(mpl108)
          		
          		After executing these commands, the header information can be accessed
          		using "mpl108.header'.  This is a Python dictionary.
          		
          		The data can be accessed using 'mpl108.dataCh1' and 'mpl108.dataCh2'.
          		These are 2-D Numpy arrays; (time, height).
            
          Tips:
               To create a pandas DataFrame of the header information, type:
                   import pandas as pd
                   header = pd.DataFrame({}, index=mpl108.header['Time'])
                   for k in mpl108.header.keys():
                       header[k] = mpl108.header[k]

"""
import os
from struct import unpack
import numpy as np
from datetime import datetime

import xarray as xr

class MPL:

    def __init__(self, filename):
        """Initialize a data object of the Micro Pulse Lidar (MPL) by reading 
        the contents of an MPL data file.  The description of the binary data 
        file is on pages 34-35 of the "Micro Pulse Lidar System, Instruction 
        Manual, MPL-4B-IDS Series" from Sigma Space Corporation.
        
        Written by Von P. Walden
        13 Oct 2014
        """

        #################################################################################
        # Open file and read critical parameters for determining header and record sizes.
        self.filename             = filename
        self.fp                   = open(filename,'rb')
        self.unitNumber           = np.uint16( unpack('<H',self.fp.read(2)))[0]
        self.fp.seek(56)
        self.numberChannels       = np.uint16( unpack('<H',self.fp.read(2)))[0]
        self.numberBins           = np.uint16( unpack('<H',self.fp.read(2)))[0]
        self.fp.seek(126)
        secondaryHdrSize          = np.uint16( unpack('<H',self.fp.read(2)))[0]
        if secondaryHdrSize==0:
            self.headerSize       = 128
        else:
            self.headerSize       = secondaryHdrSize
        self.numberBytesPerRecord = self.headerSize + (self.numberChannels*self.numberBins*4)
        
        # Determine the size of the file and the number of records.
        self.fileSize             = os.path.getsize(filename)
        self.numberRecords        = (self.fileSize / self.numberBytesPerRecord).astype('int')
        
        # Reset the file pointer to the beginning of the file.
        self.fp.seek(0)
        
        #################################################################################
        # Initialize header.
        keys = ['unitNumber',
                'version',
                'year',
                'month',
                'day',
                'hours',
                'minutes',
                'seconds',
                'Time',
                'shotsSum',
                'triggerFrequency',
                'energyMonitor',
                'detectorTemp',
                'telescopeTemp',
                'laserTemp',
                'backgroundAverage',
                'backgroundStdDev',
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
                'gpsLatitude',
                'gpsLongitude',
                'cloudBaseHeight',
                'aToDdataBadFlag',
                'dataFileVersion',
                'backgrdAverage2',
                'backgrdStdDev2',
                'mcsMode',
                'firstDataBin',
                'systemType',
                'syncPulsePerSec',
                'firstBackgrdBin',
                'secondaryHdrSize']
        self.header = dict()
        for key in keys:
            if key=='unitNumber'          or \
               key=='version'             or \
               key=='year'                or \
               key=='month'               or \
               key=='day'                 or \
               key=='hours'               or \
               key=='minutes'             or \
               key=='seconds'             or \
               key=='numberChannels'      or \
               key=='numberDataBins'      or \
               key=='scanScenarioFlag'    or \
               key=='numberBackgrdBins'   or \
               key=='firstDataBin'        or \
               key=='syncPulsePerSec'     or \
               key=='firstBackgrdBin'     or \
               key=='secondaryHdrSize':
                self.header[key] = np.array([], dtype='uint16')
            elif key=='shotsSum'          or \
                 key=='numberBins':
                self.header[key] = np.array([], dtype='uint32')
            elif key=='triggerFrequency':
                self.header[key] = np.array([], dtype='int32')
            elif key=='detectorTemp'      or \
                 key=='telescopeTemp'     or \
                 key=='laserTemp'         or \
                 key=='energyMonitor'     or \
                 key=='backgroundAverage' or \
                 key=='backgroundStdDev'  or \
                 key=='binTime'           or \
                 key=='rangeCalibration'  or \
                 key=='azimuthAngle'      or \
                 key=='elevationAngle'    or \
                 key=='compassDegrees'    or \
                 key=='polarizationV0'    or \
                 key=='polarizationV1'    or \
                 key=='gpsLatitude'       or \
                 key=='gpsLongitude'      or \
                 key=='cloudBaseHeight'   or \
                 key=='backgrdAverage2'   or \
                 key=='backgrdStdDev2':
                self.header[key] = np.array([], dtype='float32')
            elif key=='aToDdataBadFlag'   or \
                 key=='dataFileVersion'   or \
                 key=='mcsMode'           or \
                 key=='systemType':
                self.header[key] = np.array([], dtype='byte')
            elif key=='Time':
                self.header[key] = np.array([], dtype=datetime)
        
        
        #################################################################################
        # Initialize data arrays.
        self.height  = np.zeros((self.numberRecords,self.numberBins))
        self.time    = np.zeros((self.numberRecords,self.numberBins))
        self.dataCh1 = np.zeros((self.numberRecords,self.numberBins))
        self.dataCh2 = np.zeros((self.numberRecords,self.numberBins))
        
        return 
        
    def readData(self):
        """Decode the channel 1 and 2 data of an MPL data file.  
        The description of the binary data file is on pages 34-35 of the 
        "Micro Pulse Lidar System, Instruction Manual, MPL-4B-IDS Series" 
        from Sigma Space Corporation.
        
        Written by Von P. Walden
        23 Oct 2014
        
        Tips:
               To create a pandas DataFrame of the header information, type:
                   import pandas as pd
                   header = pd.DataFrame({}, index=mpl108.header['Time'])
                   for k in mpl108.header.keys():
                       header[k] = mpl108.header[k]

        Updated by Von P. Walden
        15 February 2016
        """
        
        SOL          = 299792458.      # Speed of light in m s-1
        
        for rec in range(self.numberRecords):
            
            # Read and store the header information.
                    # Decode the 128-byte header.
            header                 = self.fp.read(128)
            self.header['unitNumber']        = np.append(self.header['unitNumber']       , np.uint16( unpack('<H',header[  0:  2])[0]) )
            self.header['version']           = np.append(self.header['version']          , np.uint16( unpack('<H',header[  2:  4])[0]) )
            self.header['year']              = np.append(self.header['year']             , np.uint16( unpack('<H',header[  4:  6])[0]) )
            self.header['month']             = np.append(self.header['month']            , np.uint16( unpack('<H',header[  6:  8])[0]) )
            self.header['day']               = np.append(self.header['day']              , np.uint16( unpack('<H',header[  8: 10])[0]) )
            self.header['hours']             = np.append(self.header['hours']            , np.uint16( unpack('<H',header[ 10: 12])[0]) )
            self.header['minutes']           = np.append(self.header['minutes']          , np.uint16( unpack('<H',header[ 12: 14])[0]) )
            self.header['seconds']           = np.append(self.header['seconds']          , np.uint16( unpack('<H',header[ 14: 16])[0]) )
            self.header['Time']              = np.append(self.header['Time']             , datetime(self.header['year'][rec], self.header['month'][rec], self.header['day'][rec], self.header['hours'][rec], self.header['minutes'][rec], self.header['seconds'][rec]) )
            self.header['shotsSum']          = np.append(self.header['shotsSum']         , np.uint32( unpack('<L',header[ 16: 20])[0]) )
            self.header['triggerFrequency']  = np.append(self.header['triggerFrequency'] , np.int32(  unpack('<L',header[ 20: 24])[0]) )
            self.header['energyMonitor']     = np.append(self.header['energyMonitor']    , np.uint32( unpack('<L',header[ 24: 28])[0])/1000. )
            temp0                            = np.uint32( unpack('<L',header[ 28: 32])[0])
            temp1                            = np.uint32( unpack('<L',header[ 32: 36])[0])   # Currently not used.
            temp2                            = np.uint32( unpack('<L',header[ 36: 40])[0])
            temp3                            = np.uint32( unpack('<L',header[ 40: 44])[0])
            temp4                            = np.uint32( unpack('<L',header[ 44: 48])[0])   # Currently not used.
            self.header['detectorTemp']      = np.append(self.header['detectorTemp']     , temp0/100. )
            self.header['telescopeTemp']     = np.append(self.header['telescopeTemp']    , temp2/100. )
            self.header['laserTemp']         = np.append(self.header['laserTemp']        , temp3/100. )
            self.header['backgroundAverage'] = np.append(self.header['backgroundAverage'], np.float32(unpack('<f',header[ 48: 52])[0]) )
            self.header['backgroundStdDev']  = np.append(self.header['backgroundStdDev'] , np.float32(unpack('<f',header[ 52: 56])[0]) )
            self.header['numberChannels']    = np.append(self.header['numberChannels']   , np.uint16( unpack('<H',header[ 56: 58])[0]) )
            self.header['numberBins']        = np.append(self.header['numberBins']       , np.uint32( unpack('<L',header[ 58: 62])[0]) )
            self.header['binTime']           = np.append(self.header['binTime']          , np.float32(unpack('<f',header[ 62: 66])[0]) )
            self.header['rangeCalibration']  = np.append(self.header['rangeCalibration'] , np.float32(unpack('<f',header[ 66: 70])[0]) )
            self.header['numberDataBins']    = np.append(self.header['numberDataBins']   , np.uint16( unpack('<H',header[ 70: 72])[0]) )
            self.header['scanScenarioFlag']  = np.append(self.header['scanScenarioFlag'] , np.uint16( unpack('<H',header[ 72: 74])[0]) )
            self.header['numberBackgrdBins'] = np.append(self.header['numberBackgrdBins'], np.uint16( unpack('<H',header[ 74: 76])[0]) )
            self.header['azimuthAngle']      = np.append(self.header['azimuthAngle']     , np.float32(unpack('<f',header[ 76: 80])[0]) )
            self.header['elevationAngle']    = np.append(self.header['elevationAngle']   , np.float32(unpack('<f',header[ 80: 84])[0]) )
            self.header['compassDegrees']    = np.append(self.header['compassDegrees']   , np.float32(unpack('<f',header[ 84: 88])[0]) )
            self.header['polarizationV0']    = np.append(self.header['polarizationV0']   , np.float32(unpack('<f',header[ 88: 92])[0]) )
            self.header['polarizationV1']    = np.append(self.header['polarizationV1']   , np.float32(unpack('<f',header[ 92: 96])[0]) )
            self.header['gpsLatitude']       = np.append(self.header['gpsLatitude']      , np.float32(unpack('<f',header[ 96:100])[0]) )
            self.header['gpsLongitude']      = np.append(self.header['gpsLongitude']     , np.float32(unpack('<f',header[100:104])[0]) )
            self.header['cloudBaseHeight']   = np.append(self.header['cloudBaseHeight']  , np.float32(unpack('<f',header[104:108])[0]) )
            self.header['aToDdataBadFlag']   = np.append(self.header['aToDdataBadFlag']  , np.byte(   unpack('<b',header[108:109])[0]) )
            self.header['dataFileVersion']   = np.append(self.header['dataFileVersion']  , np.byte(   unpack('<b',header[109:110])[0]) )
            self.header['backgrdAverage2']   = np.append(self.header['backgrdAverage2']  , np.float32(unpack('<f',header[110:114])[0]) )
            self.header['backgrdStdDev2']    = np.append(self.header['backgrdStdDev2']   , np.float32(unpack('<f',header[114:118])[0]) )
            self.header['mcsMode']           = np.append(self.header['mcsMode']          , np.byte(   unpack('<b',header[118:119])[0]) )
            self.header['firstDataBin']      = np.append(self.header['firstDataBin']     , np.uint16( unpack('<H',header[119:121])[0]) )
            self.header['systemType']        = np.append(self.header['systemType']       , np.byte(   unpack('<b',header[121:122])[0]) )
            self.header['syncPulsePerSec']   = np.append(self.header['syncPulsePerSec']  , np.uint16( unpack('<H',header[122:124])[0]) )
            self.header['firstBackgrdBin']   = np.append(self.header['firstBackgrdBin']  , np.uint16( unpack('<H',header[124:126])[0]) )
            self.header['secondaryHdrSize']  = np.append(self.header['secondaryHdrSize'] , np.uint16( unpack('<H',header[126:128])[0]) )
            
            # Skips the bytes at the end of the header; extra space for secondary header.
            # This is only true for MPL 108 (ARM MPL for N-ICE 2015), not MPL 107, which
            #     is being used at Summit Station as part of the ICECAPS experiment.
            if self.header['unitNumber'][rec]==108:
                self.fp.read(self.header['secondaryHdrSize'][rec]-128)
            
            # Read and store the data records.
            self.dataCh1[rec,:] = np.fromstring(self.fp.read(4*self.numberBins), dtype='<f4')
            self.dataCh2[rec,:] = np.fromstring(self.fp.read(4*self.numberBins), dtype='<f4')
            
        
        # Calculates the height and time vectors.
        rng         = 0.5 * SOL * self.header['binTime'][0] * 0.001	# Range gate altitude in km
        self.height = (range(self.header['numberBins'][0]) - self.header['firstDataBin'][0]) * rng
        self.time   = self.header['Time']
        self.hours  = np.array([])
        for t in self.time:
            self.hours = np.append(self.hours, (t - datetime(self.header['year'][0],self.header['month'][0],self.header['day'][0])).total_seconds()/(3600.))
            
        self.fp.close()
        
        return
    
    def to_xarray(self):
        header = {k: ('time', self.header[k]) for k in self.header.keys()}
        data   = {'backscatter_1': (('time', 'height'), self.dataCh1[:,0:1200]), 'backscatter_2': (('time', 'height'), self.dataCh2[:,0:1200])}
        ds = xr.Dataset(
            data_vars= header | data, 
            coords={'time': self.time, 'height': self.height[0:1200]}
        )

        return ds

    def plotData(self):
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        
        plt.figure()
        plt.subplot(311)
        plt.imshow(np.flipud(np.transpose(self.dataCh1[:,200:])),extent=[min(self.hours),max(self.hours),min(self.height[200:]),max(self.height[200:])],aspect='auto',cmap='viridis',norm=LogNorm())
        plt.colorbar()
        plt.title('ARM MPL field data (Channel 1), ' + self.time[0].strftime('%d %b %Y') + ', RV Lance')
        
        plt.subplot(312)
        plt.imshow(np.flipud(np.transpose(self.dataCh2[:,200:])),extent=[min(self.hours),max(self.hours),min(self.height[200:]),max(self.height[200:])],aspect='auto',cmap='viridis',norm=LogNorm())
        plt.colorbar()
        plt.ylabel('Height (km)')
        plt.title('ARM MPL field data (Channel 2), ' + self.time[0].strftime('%d %b %Y') + ', RV Lance')
        
        plt.subplot(313)
        plt.imshow(np.flipud(np.transpose(self.dataCh1[:,200:]/self.dataCh2[:,200:])),extent=[min(self.hours),max(self.hours),min(self.height[200:]),max(self.height[200:])],aspect='auto',cmap='viridis',norm=LogNorm())
        plt.colorbar()
        plt.xlabel('Time of Day (UTC hours)')
        plt.title('ARM MPL field data (depol ratio), ' + self.time[0].strftime('%d %b %Y') + ', RV Lance')
                
        plt.show()
        return

    def getDate(self):

        YY = int(self.filename[-16:-12])
        mm = int(self.filename[-12:-10])
        dd = int(self.filename[-10:-8])

        self.date = datetime(YY,mm,dd)
        return
    
    def plotNICEdata(self):
        """
        One can simply initialize mpl.MPL with any file from a day and plotNICEdata will
        plot ALL the data from that day.
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from datetime import datetime
        from glob import glob

        self.getDate()
        day = self.date
        fns = np.sort(glob(self.filename[:-8] + '*.mpl'))
        if len(fns)==0:
            print ('Missing day: %s') % day.strftime('%Y-%m-%d')
        else:
            print(fns[0])
            mpl108 = mpl.MPL(fns[0])
            mpl108.readData()
            ch1    = mpl108.dataCh1
            ch2    = mpl108.dataCh2
            Time   = np.array([])
            for t in mpl108.header['Time']:
                Time = np.append(Time, (t-datetime(2015,4,24)).seconds/3600.)
            for fn in fns[1:]:
                print(fn)
                mpl108 = mpl.MPL(fn)
                mpl108.readData()
                ch1    = np.concatenate((ch1, mpl108.dataCh1), axis=0)
                ch2    = np.concatenate((ch2, mpl108.dataCh2), axis=0)
                for t in mpl108.header['Time']:
                    Time = np.append(Time, (t-datetime(2015,4,24)).seconds/3600.)
        
        # Plot the data for the desired day.
        plt.figure()
        plt.subplot(311)
        plt.xticks(np.array([]))
        plt.imshow(np.flipud(np.transpose(ch1[:,200:868])),extent=[min(Time),max(Time),mpl108.height[200],mpl108.height[868]],aspect='auto',cmap='viridis',norm=LogNorm())
        plt.axis([0.,24.,0.,10.])
        plt.colorbar()
        plt.title('CHANNEL 1: ARM MPL field data, ' + day.strftime('%d %b %Y') + ', RV Lance')
        
        plt.subplot(312)
        plt.xticks(np.array([]))
        plt.imshow(np.flipud(np.transpose(ch2[:,200:868])),extent=[min(Time),max(Time),mpl108.height[200],mpl108.height[868]],aspect='auto',cmap='viridis',norm=LogNorm())
        plt.axis([0.,24.,0.,10.])
        plt.colorbar()
        plt.ylabel('Height (km)')
        plt.title('CHANNEL 2')
        
        plt.subplot(313)
        plt.xticks(np.arange(0.,25.,6.))
        plt.imshow(np.flipud(np.transpose(ch1[:,200:868]/ch2[:,200:868])),extent=[min(Time),max(Time),mpl108.height[200],mpl108.height[868]],aspect='auto',cmap='viridis',norm=LogNorm())
        plt.axis([0.,24.,0.,10.])
        plt.colorbar()
        plt.xlabel('Time of Day (UTC hours)')
        plt.title('DEPOLARIZATION RATIO')
                
        plt.savefig('/Volumes/N-ICE/WP2/MPL/DATA/images/mpl_'+ day.strftime('%Y%m%d') + '.png')
        plt.close()
        return
    
    def plotTemperatureData(self):
        import matplotlib.pyplot as plt
        
        plt.figure()
        plt.plot(self.time,self.header['telescopeTemp'],self.time, self.header['laserTemp'],self.time,self.header['detectorTemp'])
        plt.xlabel('Time (UTC)')
        plt.ylabel('Temperature (C)')
        plt.title('ARM MPL field data (Instrument Temperatures), ' + self.time[0].strftime('%d %b %Y') + ', RV Lance')
        plt.legend(('telescope','laser','detector'),loc='best')
        
        plt.show()
        return
    
    def plotEnergyMonitorData(self):
        import matplotlib.pyplot as plt
        
        plt.figure()
        plt.plot(self.time,self.header['energyMonitor'])
        plt.xlabel('Time (UTC)')
        plt.ylabel('Laser Energy (mJ)')
        plt.title('ARM MPL field data (Energy Monitor), ' + self.time[0].strftime('%d %b %Y') + ', RV Lance')
        
        plt.show()
        return
    
