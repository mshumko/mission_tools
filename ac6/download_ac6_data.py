# This code downloads AC6 data for the days specified
import requests
import urllib
import tarfile
import io
from datetime import datetime, timedelta
# import itertools

def download_ac6_day(sc_id, date, saveDir=None):
    """
    This function downloads AC-6 data from the RBSP science gateway and 
    saves it to a saveDir. The tgz file is uncompressed.
    """
    print( date.year)
    url = ('http://rbspgway.jhuapl.edu/share/ac6/data/'
        'AC6-{}/{}/AC6-{}_{}_V03.tgz'.format(
            sc_id.upper(), date.year, sc_id.upper(), date.date().isoformat().replace('-', '')))
    print(url)
    #request = requests.get(url) # Get tgz file
    request = urllib.request.urlopen(url)
    #print(request)
    #tar = tarfile.open(io.BytesIO(request.read()), 'r')
    #print(tar)
    return request

if __name__ == '__main__':
    import argparse
    import itertools
    parser = argparse.ArgumentParser(description=('This script downloads the '
            'AC6 data from a website defined in the code.'))
    parser.add_argument('-sc', '--sc_id', help=('Define a spaceraft id. '
                                            'Use either "a" or "b". If '
                                            'none, will download data '
                                            'from both AC-6 spacecraft.'),
                                            type=str, nargs=1)
    parser.add_argument('date', default=[], nargs='*', type=int,
                        help=('This is an optional date positional argument.'
                                ' The format is YYYY MM DD.'))
    parser.add_argument('-dr', '--dateRange', default=[], nargs='*', type=int,
                            help=('This optional argument passes in a start'
                                    ' and end dates to define a date range for '
                                    'which the script will download EMFISIS data.'
                                    ' The format is YYYY MM DD and the input is '
                                    'size 6'))
    args = parser.parse_args()
    if args.date != []:
        if len(args.date) != 3:
            raise ValueError('Date not formated correctly! Check manual.')
        startDate = datetime(*args.date)
        endDate = datetime(*args.date)
    if args.dateRange != []:
        if len(args.dateRange) != 6:
            raise ValueError('Date range not formated correctly! Check manual.')
        startDate = datetime(*args.dateRange[0:3])
        endDate = datetime(*args.dateRange[3:])
    if args.sc_id is None:
        args.sc_id = ['A', 'B']

    # Loop over the dates:
    dates = [startDate + timedelta(days=d) for d in 
        range(1+(endDate - startDate).days)]

    for (sc_id, day) in itertools.product(args.sc_id, dates):
        r  =download_ac6_day(sc_id, day)
