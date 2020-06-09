# This program loads the HILT data and parses it into a nice format
import argparse
import pandas as pd
import pathlib
from datetime import datetime, date
import zipfile
import numpy as np

hilt_dir = '/home/mike/research/sampex/'

class Load_HILT:
    def __init__(self, load_date, zipped=True, extract=False, 
                time_index=True):
        """
        Load the HILT data given a date. If zipped is True, this class will
        look for txt.zip file with the date and open it (without extracting).
        If you want to extract the file as well, set extract=True.
        time_index=True sets the time index of self.hilt to datetime objects
        otherwise the index is just an enumerated list.
        """
        self.load_date = load_date
        if zipped:
            extention='.txt.zip'
        else:
            extention='.txt'
        
        # Figure out how to calculate the day of year (DOY)
        if isinstance(self.load_date, pd.Timestamp):
            doy = str(self.load_date.dayofyear).zfill(3)
        elif isinstance(self.load_date, (datetime, date) ):
            doy = str(self.load_date.timetuple().tm_yday).zfill(3)

        # Get the filename and search for it. If multiple or no
        # unique files are found this will raise an assertion error.
        file_name = f'hhrr{self.load_date.year}{doy}{extention}'
        matched_files = list(pathlib.Path(hilt_dir).rglob(file_name))
        assert len(matched_files) == 1, (f'0 or >1 matched HILT files found.'
                                        f'\n{file_name}'
                                        f'\nmatched_files={matched_files}')
        self.file_path = matched_files[0]

        # Load the zipped data and extract if it is set to true
        if zipped:
            self.read_zip(self.file_path, extract=extract)
        else:
            self.read_csv(self.file_path)

        # Parse the seconds of day time column to datetime objects
        self.parse_time(time_index=time_index)
        return
        
    def read_zip(self, zip_path, extract=False):
        """
        Open the zip file and load in the csv file. If extract=False than the file
        will only be opened and not extracted to a text file in the 
        sampex/data/hilt directory. 
        """
        txt_name = zip_path.stem # Remove the .zip from the zip_path
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            if extract:
                zip_ref.extractall(zip_path.parent)
                #self.hilt = pd.read_csv(zip_path.parent / txt_name)
                self.read_csv(zip_path.parent / txt_name)
            else:
                with zip_ref.open(txt_name) as f:
                    # self.hilt = pd.read_csv(f, sep=' ')
                    self.read_csv(f)
        return
    
    def read_csv(self, path):
        """
        Reads in the CSV file given either the filename or the 
        zip file reference
        """
        print(f'Loading SAMPEX on {self.load_date.date()} from {path.name}')
        self.hilt = pd.read_csv(path, sep=' ')
        return

    def parse_time(self, time_index=True):
        """ 
        Parse the seconds of day column to a datetime column. 
        If time_index=True, the time column will become the index.
        """
        day_seconds_obj = pd.to_timedelta(self.hilt['Time'], unit='s')
        self.hilt['Time'] = pd.Timestamp(self.load_date) + day_seconds_obj
        if time_index:
            self.hilt.index = self.hilt['Time']
            del(self.hilt['Time'])
        return

class Load_Attitude:
    def __init__(self, date):

        return

if __name__ == '__main__':
    l = Load_HILT(datetime(2000, 4, 4))

    #### Use argparse if running in interactive mode ###
    # parser = argparse.ArgumentParser(description=('This script plots the '
    #     'SAMPEX HILT data.'))
    # parser.add_argument('date', nargs=3, type=int,
    #     help=('This is the date to plot formatted as YYYY MM DD')) 
    # parser.add_argument('-d', '--dtype', type=str, default='10Hz',
    #     help=('AC6 data type to plot (10Hz or survey)'))  
    # parser.add_argument('-p', '--plot', type=bool, default=True,
    #     help=('Plot AC6 data'))
    # args = parser.parse_args()

    # date = datetime(*args.date)
    # import time
    # t = time.time()
    # data = read_ac_data_wrapper(args.sc_id, date, dType=args.dtype, 
    #         tRange=None)

    # if args.plot:
    #     p = Plot_AC6(data, args.sc_id, args.dtype)
    #     p.plot_data()