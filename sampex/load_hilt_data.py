# This program loads the HILT data and parses it into a nice format
import argparse
import pandas as pd
import pathlib
from datetime import datetime, date
import zipfile
import numpy as np
import re

hilt_dir = '/home/mike/research/sampex/data/hilt'
attitude_dir = '/home/mike/research/sampex/data/attitude'

class Load_SAMPEX_HILT:
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
        print(f'Loading SAMPEX HILT data from {self.load_date.date()} from {path.name}')
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

    def resolve_counts_state4(self):
        """ 
        This function resolves the HILT counts to 20 ms resolution assuming 
        the data is in state4. The counts represent the sum from the 4 SSDs.
        Data saved in self.hilt_resolved
        """ 
        resolution_ms = 20E-3
        # Resolve the counts using numpy (most efficient way with 
        # static memory allocation)
        self.counts = np.nan*np.zeros(5*self.hilt.shape[0], dtype=int)
        for i in [0, 1, 2, 3]:
            self.counts[i::5] = self.hilt[f'Rate{i+1}']
        # This line is different because rate5 is 100 ms SSD4 data.
        self.counts[4::5] = self.hilt['Rate6'] 

        # Resolve the time array.
        self.times = np.nan*np.zeros(5*self.hilt.shape[0], dtype=object)
        for i in [0, 1, 2, 3, 4]:
            self.times[i::5] = self.hilt.index + pd.to_timedelta(resolution_ms*i, unit='s')

        self.hilt_resolved = pd.DataFrame(data={'counts':self.counts}, index=self.times)
        return self.counts, self.times


class Load_SAMPEX_Attitude:
    def __init__(self, load_date):
        """ 
        This class loads the appropriate SAMEX attitude file, 
        parses the complex header and converts the time 
        columns into datetime objects
        """
        self.load_date = load_date

        # Figure out how to calculate the day of year (DOY)
        if isinstance(self.load_date, pd.Timestamp):
            self.doy = int(self.load_date.dayofyear)
        elif isinstance(self.load_date, (datetime, date) ):
            self.doy = int(self.load_date.timetuple().tm_yday)

        # Find the appropriate attitude file.
        self.find_matching_attitude_file()

        # Load the data into a dataframe
        self.load_attitude()
        return

    def find_matching_attitude_file(self):
        """ 
        Uses pathlib.rglob to find the attitude file that contains 
        the DOY from self.load_date
        """
        attitude_files = list(pathlib.Path(attitude_dir).rglob('PSSet_6sec_*_*.txt'))
        start_end_dates = [re.findall(r'\d+', str(f))[1:] for f in attitude_files]
        
        current_date_int = int(self.load_date.year*1000 + self.doy)
        self.attitude_file = None

        for f, (start_date, end_date) in zip(attitude_files, start_end_dates):
            if (int(start_date) <= current_date_int) and (int(end_date) >= current_date_int):
                self.attitude_file = f
        if self.attitude_file is None:
            raise ValueError(f'A matched file not found for year='
                             f'{self.load_date.year}, doy={self.doy}')
        return self.attitude_file

    def load_attitude(self, columns='default', remove_old_time_cols=True):
        """ 
        Loads the attitude file. Only columns specified in the columns arg are 
        loaded to conserve memory. The year, day_of_year, and sec_of_day columns
        are used to construct a list of datetime objects that are assigned to
        self.attitude index. 

        If remove_old_time_cols is True, the year, DOY, and second columns are 
        delited to conserve memory.
        """
        print(f'Loading SAMPEX attitude data from {self.load_date} from'
            f' {self.attitude_file.name}')
        # A default set of hard-coded list of columns to load
        if columns=='default':
            columns = {
                0:'Year', 1:'Day-of-year', 2:'Sec_of_day', 6:'GEO_Radius',
                7:'GEO_Long', 8:'GEO_Lat', 9:'Altitude', 20:'L_Shell',
                22:'MLT', 42:'Mirror_Alt', 68:'Pitch'
            }
        # Open the attitude file stream
        with open(self.attitude_file) as f:
            # Skip the long header until the "BEGIN DATA" line 
            self._skip_header(f)
            # Save the rest to a file using columns specified by the columns.keys() with the 
            # columns values for the column names.
            self.attitude = pd.read_csv(f, sep=' ', names=[columns[key] for key in columns.keys()], 
                                        usecols=columns.keys())
        # Parse the dates by first making YYYY-DOY strings.
        year_doy = [f'{year}-{doy}' for year, doy in 
                    self.attitude[['Year', 'Day-of-year']].values]
        attitude_dates=pd.to_datetime(year_doy, format='%Y-%j')
        # Now add the seconds of day to complete the date and time.
        self.attitude.index = attitude_dates + pd.to_timedelta(self.attitude['Sec_of_day'], unit='s')
        # Optionally remove duplicate columns to conserve memory.
        if remove_old_time_cols:
            self.attitude.drop(['Year', 'Day-of-year', 'Sec_of_day'], axis=1, inplace=True)
        return

    def _skip_header(self, f):
        """ 
        Read in the "f" attitude file stream line by line until the 
        "BEGIN DATA" line is reached. Then return Returns a list of column
        names from the parsed header.
        """
        for line in f:
            if "BEGIN DATA" in line:
                return f 
        return None


class Load_SAMPEX_HILT_ATTITUDE:
    def __init__(self, date):
        """
        This is a child class of Load_SAMPEX_HILT and Load_SAMPEX_Attitude
        that loads the SAMPEX HILT and attitude data and merges the two
        datasets using the merge_asof pandas function.
        """
        raise NotImplementedError
        return

    def merge_data(self):
        """ 
        This method uses pd.merge_asof to merge the SAMPEX attitude data with
        the HILT data with a maximum threshold of 6 seconds.
        """ 
        raise NotImplementedError
        return

if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt

    start_time = time.time()

    l = Load_SAMPEX_HILT(datetime(2000, 4, 4))
    l.resolve_counts_state4()
    a = Load_SAMPEX_Attitude(datetime(2000, 4, 4))

    print(f'Run time = {time.time()-start_time} s')

    plt.plot(l.hilt_resolved.index, l.hilt_resolved.counts)
    plt.show()


    ### Copied AC6 code to make a command-line interface to plot the daily SAMPEX data 
    ### once Load_SAMPEX_HILT_ATTITUDE is written.
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