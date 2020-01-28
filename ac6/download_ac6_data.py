import os
import urllib3
import tarfile
import io
import glob
import argparse
from datetime import datetime, date, timedelta

def download_ac6_data(sc_id, data_date, save_dir, verbose=True):
    """
    This function downloads AC6 data.

    Parameters
    ----------
    sc_id : str
        Spacecraft id, either A or B. It is case insensitive.
    data_date : str, datetime, or date object
        The date you need the path for. The date can a string 
        formatted as, for example, 20161010 for 10 Oct 2016 
        date. The date can also be a datetime or date object.
    save_dir : str
        The path where to save the AC6 csv files.
    verbose : bool, optional
        If true will print a status message if data was or was not
        downloaded for that day.

    Returns
    -------
    None
    """
    file_path = get_ac6_path(sc_id, date)

    http = urllib3.PoolManager()
    f = http.request('GET', file_path)
    tarball = io.BytesIO(f.data)

    try:
        # If tarball file found, extract to save_dir and exit.
        tar = tarfile.open('r', fileobj=tarball)
        tar.extractall(save_dir)
        # Rename the files from the YYYYDDMM convention to YYYYMMDD
        rename_files(sc_id, data_date, save_dir) 
        if verbose: print(f'Data for AC6-{sc_id.upper()} on {data_date} downloaded')
    except tarfile.ReadError:
        # If file not found, print error if verbose and exit.
        if verbose: print(f'Data for AC6-{sc_id.upper()} on {data_date} could not '
                          'be downloaded')
    return

def get_ac6_path(sc_id, data_date, 
                base_url='http://rbspgway.jhuapl.edu/share/ac6/data/'):
    """
    This function returns a url path for the AC6 data given 
    spacecraft id and date. The date can be either a string
    formatted as 20161010 for 10 Oct 2016 date.

    Parameters
    ----------
    sc_id : str
        Spacecraft id, either A or B. It is case insensitive.
    data_date : str, datetime, or date object
        The date you need the path for. The date can a string 
        formatted as, for example, 20161010 for 10 Oct 2016 
        date. The date can also be a datetime or date object.

    Returns
    -------
    file_path : str 
        A complete file path (url) that specifies where the AC6 
        data lives on the RBSP science gateway for spacecraft 
        sc_id on date data_date.
    """
    if isinstance(data_date, str):
        # if data_date formatted as 20161010 for 10 Oct 2016 date.
        year = data_date[:4]
        file_dir = f'AC6-{sc_id.upper()}/{year}/'
        file_name = f'AC6-{sc_id.upper()}_{data_date}_V03.tgz'
    elif isinstance(data_date, datetime) or isinstance(data_date, date):
        # Otherwise if the data_date is a datetime or date object.
        file_dir = f'AC6-{sc_id.upper()}/{data_date.year}/'
        file_name = f'AC6-{sc_id.upper()}_{data_date.strftime("%Y%d%m")}_V03.tgz'
    else:
        raise ValueError('The data_data argument must be a string,'
                        ' datetime, or date object.')
    file_path = os.path.join(base_url, file_dir, file_name)
    return file_path

def rename_files(sc_id, data_date, save_path):
    """
    This function rename the downloaded AC6 files from the 
    YYYYDDMM convention to YYYYMMDD

    Parameters
    ----------
    sc_id : str
        Spacecraft id, either A or B. It is case insensitive.
    data_date : str, datetime, or date object
        The date you need the path for. The date can a string 
        formatted as, for example, 20161010 for 10 Oct 2016 
        date. The date can also be a datetime or date object.
    save_dir : str
        The path where to save the AC6 csv files.

    Returns
    -------
    None
    """
    old_file_search_string = (f'AC6-{sc_id.upper()}_'
                            f'{data_date.strftime("%Y%d%m")}_*.csv')
    old_files = glob.glob(os.path.join(save_path, old_file_search_string))

    for old_file_path in old_files:
        split_path = old_file_path.split('_')
        split_path[-4] = data_date.strftime('%Y%m%d')
        new_file_path = '_'.join(split_path)
        os.rename(old_file_path, new_file_path)
    return


if __name__ == '__main__':
    """
    This script wrapps the download_ac6_data function and uses command line
    arguments to download AC6 data from multiple days.

    For example, if you want data from both AC6 units from 1 October 2016 to
    31 October 2016 call:

    python3 download_ac6_data.py A B 2016 10 1 2016 10 31
    """

    parser = argparse.ArgumentParser(
                    description=('This script downloads AC6 data from '
                                'the RBSP science gateway.'))
    parser.add_argument('sc_id', type=str, nargs='+',
        help=('Spacecraft to download data from. '
        'Can be "A", "B", or "A B".'))
    parser.add_argument('dates', nargs=6, type=int,
        help=('Start and end dates to download. Must be in the following format '
        'YYYY MM DD YYY MM DD.')) 
    parser.add_argument('-sd', '--save_dir', type=str, 
        default='/home/mike/Downloads/ac6', 
        help=('Directory where to save the data to.'))
    args = parser.parse_args()

    # Figure out which days to download the data from.
    start_date = datetime(*args.dates[:3])
    end_date = datetime(*args.dates[3:])
    dates = [start_date + timedelta(days=i) 
            for i in range((end_date-start_date).days+1)]

    for date in dates:
        for sc_id in args.sc_id:
            download_ac6_data(sc_id, date, args.save_dir)


