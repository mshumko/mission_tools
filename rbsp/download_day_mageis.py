# Download MagEIS public Release 3 data
import urllib.request, os
from datetime import datetime
from datetime import timedelta
import itertools
import re

def downloadMagEIS(rb_id, dates, fPath, dataLevel=2):
    """
    This function takes in an rbsp id, a date or array of dates, save directory
    and a data level and downloads data from the LANL website.
    """
    # Allow the dates input to be a single date or an array of dates.
    if not hasattr(dates, '__len__'):
        dates = [dates]
        
    for date in dates:
        print('Downloading RBSP-{} MagEIS data from {}'.format(rb_id.upper(), date.date()))
        if dataLevel == 2:
            url = ('https://rbsp-ect.lanl.gov/data_pub/rbsp{}/mageis/level2/sectors/'
                '{}/'.format(rb_id.lower(), date.year))
        else:
            url = ('https://rbsp-ect.lanl.gov/data_pub/rbsp{}/mageis/level3/'
                'pitchangle/{}/'.format(rb_id.lower(), date.year))
        
        # Find the filename(s) The version numbers are different, so need to use
        # regular expressions
        joinDate = ''.join(date.date().isoformat().split('-'))
        response = urllib.request.urlopen(url)
        data = response.read().decode()  # Get HTML data to find filename
        fName = findFilename(rb_id, joinDate, data)      
        if fName is None:  
            continue
            
        # Download data
        urllib.request.urlretrieve(url + fName, os.path.join(fPath, fName))
    return 
    
def findFilename(rb_id, joinDate, data):
    """
    This function uses regular expressions to find the correct
    filename to download
    """
    regex = re.compile(r'rbsp{}_rel\d\d_ect-mageis-L\d_{}_v'
        '\d.\d.\d.cdf'.format(rb_id.lower(), joinDate))
    mo = regex.search(data)
    try:
        return mo.group()
    except AttributeError as err:
        print(err, 'Probably no data was found on this date.')
        return None
    
if __name__ == '__main__':
    import argparse
    import itertools
    parser = argparse.ArgumentParser(description=('This script downloads the '
            'MagEIS data from a website defined in the code.'))
    parser.add_argument('-sc', '--sc_id', help=('Define a spaceraft id. '
                                            'Use either "a" or "b". If '
                                            'none, will download data '
                                            'from both RBSP spacecraft.'),
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


    # START_DATE = datetime(2017, 11, 18)
    # END_DATE = datetime(2017, 11, 19)
    dates = [startDate + timedelta(days=d) for d in 
        range(1+(endDate - startDate).days)]
    #for (rb_id, day) in itertools.product(['a', 'b'], dates):
    for rb_id in args.sc_id:
        fPath = '/home/mike/research/rbsp/data/mageis/rbsp{}'.format(rb_id.lower())
        downloadMagEIS(rb_id, dates, fPath, dataLevel=3)
    

