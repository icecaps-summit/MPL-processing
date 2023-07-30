'''Author: Andrew Martin
Date created: 15/2/23

Python script for creating ingested mpl files for a given month.
'''

import datetime
import xarray as xr
import numpy as np
import os

from . import load
from . import raw_to_ingested
from . import afterpulse
from . import calibrate_ingested

def create_ingested(date,dir_target,dir_mpl, overwrite=False, afterpulse=None, overlap=None, sources={}):
    '''For a given date, create an ingested file and save it to the target directory.
    
    INPUTS:
        date : datetime.date
            The date for which data should be taken from.
            
        dir_target : string
            The directory into which the file should be saved.

        dir_mpl : string
            The directory from which the .mpl.gz files will be extracted.
            
        overwrite : boolean
            If true, pre-existing ingested files will be overwritten.

        afterpulse : xr.Dataset
            xarray dataset containing the information from the afterpulse file generated by load_afterpusle

        fname_overlap : xr.Dataset | 2xk np.ndarray
            Data for the overlap function that can be utilised in calibrate_ingested.
    '''
    # check to see if the ingested file already exists, and if it can be overwritten.
    #save_fname = f'smtmplpolX1.a1.{date.year:04}{date.month:02}{date.day:02}.000000.cdf'
    save_fname = f'mpl_ingested_{date.year:04}{date.month:02}{date.day:02}.nc'
    if not overwrite:
        if os.path.isfile(os.path.join(dir_target,save_fname)):
            print(f'{save_fname} already exists.')
            return

    # create the filename format from the date
    fname_glob = f'{date.year:04}{date.month:02}{date.day:02}*00.mpl.gz'
    print(fname_glob)
    ds = load.load_fromglob(fname_glob, dir_mpl)

    # apply the raw_to_ingested algorithm on the already-loaded ds
    ds = raw_to_ingested.raw_to_ingested(None, None, data_loaded=ds)

    # add calibrated variables to the ingested format
    ds = calibrate_ingested.calibrate_ingested(ds, afterpulse=afterpulse, overlap=overlap, sources=sources)

    # now save the dataset as a netcdf file
    ds.to_netcdf(os.path.join(dir_target, save_fname))
    return



def create_ingested_month(year, month, dir_target, dir_mpl, overwrite=False, fname_afterpulse=None, fname_overlap=None):
    '''Function to call create_ingested() for every day in a month.
    
    INPUTS:
        year : int
            Year for data to be converted.
            
        month : int
            Month for data to be converted.

        dir_target : string
            Directory for the ingested files to be saved to.

        dir_mpl : string
            The directory from which the .mpl.gz files will be extracted.

        overwrite : boolean
            If true, pre-existing ingested files can be overwritten.

        fname_afterpulse : string
            Full filename for the afterpulse file.

        fname_overlap : string
            Full filename for the file containing the overlap function data.
    '''
    overlap, afterpulse, sources = load_o_a_s(fname_overlap, fname_afterpulse)

    delta_day = datetime.timedelta(days=1)
    date0 = datetime.date(year=year, month=month, day=1)
    while date0.month == month:
        create_ingested(date=date0, dir_target=dir_target, dir_mpl=dir_mpl, overwrite=overwrite, afterpulse=afterpulse, overlap=overlap, sources=sources)
        date0 = date0 + delta_day # increment the date



def load_o_a_s(fname_overlap,fname_afterpulse):
    '''Function to laod in the overlap and afterpulse data
    
    INPUTS:
        fname_overlap : string, None
            Full filename for the overlap data. If None, will default to Turner values.
            
        fname_afterpulse : string, None
            Full filename for the afterpulse data.
            
    OUTPUTS:
        overlap : xr.Dataset | (2,k) np.ndarray
            Object containing the overlap function data for use in calibrate_ingested.
            
        afterpulse : xr.Dataset
            xarray dataset containing the afterpulse data for use in calibrate_ingested.
            
        sources : dictionary
            Dictionary for containing the data sources for the overlap and afterpulse data.
    '''
    sources = {}
    # load in the overlap function
    if fname_overlap is not None:
        sources['overlap'] = fname_overlap
        raise NotImplementedError('loading of overlap function not supported yet.')
    else:
        sources['overlap'] = 'David Turner values from 31/01/2011.'
        ocorr = np.array([0.00530759, 0.0583835, 0.110524, 0.174668, 0.246036, 0.333440, 0.421466,0.510560, 0.599191, 0.676644, 0.744512, 0.808004, 0.848976,0.890453, 0.959738, 0.975302, 1.0, 1.0])
        oht = np.array([0.0149896, 0.164886, 0.314782, 0.464678, 0.614575, 0.764471, 0.914367,1.06426, 1.21416, 1.36406, 1.51395, 1.66385, 1.81374, 1.96364,2.11354, 2.26343, 2.5, 20]) * 1e3
        overlap = np.vstack((oht, ocorr))

    # load in the afterpulse file
    afterpulse = None
    if fname_afterpulse is not None:
        sources['afterpulse'] = fname_afterpulse
        afterpulse = afterpulse.load_afterpulse(fname_afterpulse)

    return overlap, afterpulse, sources



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to convert archived .mpl.gz files to ingested .cdf files for a given month.')

    parser.add_argument('-y', '--year', required=True, type=int, help='The year for the month being converted.')
    parser.add_argument('-m', '--month', required=True, type=int, help='Month for the data being converted, as an integer.')
    parser.add_argument('-t', '--targetdir', default='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/mpl/leeds_ingested', help='The directory that the ingested files will be saved to. Defaults to /gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/mpl/leeds_ingested')
    parser.add_argument('-d', '--datadir', default='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/mpl/raw', help='The directory from which the raw .mpl.gz data will be extracted. Defaults to /gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/mpl/raw')
    parser.add_argument('-o', '--overwrite', action='store_true', help='Optional, Overwrite existing ingested files at targetdir.')

    parser.add_argument('-A', '--afterpulse', help='Optional, Full filename for the afterpulse file.')
    parser.add_argument('-O', '--overlap', help='Optional, Full filename for the overlap function file.')

    # an optional argument, if day is passed in then we just do a single day
    parser.add_argument('--day', type=int, help='Optional, specifies a particular day for which the ingestion should be done.')

    args = parser.parse_args()

    year = args.year
    month = args.month
    dir_target = args.targetdir
    dir_mpl = args.datadir
    overwrite = args.overwrite

    fname_afterpulse = args.afterpulse
    fname_overlap = args.overlap

    day = args.day
    if day is not None:
        date0 = datetime.date(year=year, month=month, day=day)

        overlap, afterpulse, sources = load_o_a_s(fname_overlap=fname_overlap, fname_afterpulse=fname_afterpulse)
        create_ingested(date=date0, dir_target=dir_target, dir_mpl=dir_mpl, overwrite=overwrite, afterpulse=afterpulse, overlap=overlap, sources=sources)
    else:
        create_ingested_month(year=year, month=month, dir_target=dir_target, dir_mpl=dir_mpl, overwrite=overwrite, fname_afterpulse=fname_afterpulse, fname_overlap=fname_overlap)
