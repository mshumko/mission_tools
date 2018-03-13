import urllib.request, os
from datetime import datetime
from datetime import timedelta
import re
#url = 'https://rbsp-ect.lanl.gov/data_pub/rbspa_bak/MagEphem/def/2016/rbspa_def_MagEphem_T89D_20160105_v1.0.0.txt

def findRBSPMagEphemUrl(sc_id, date, **kwargs):
    """
    sc_id is either 'A' or 'B'
    """
    url = kwargs.get('url', None)
    Bmodel = kwargs.get('Bmodel', 'T89D')
    fType = kwargs.get('fType', 'txt')
    splitDate = date.date().isoformat().split('-')

    # Set up the regular expression to find the file.
    reString = r'rbsp{}.*{}.*{}.*{}'.format(sc_id.lower(), Bmodel, ''.join(splitDate), fType)
    magEphemRegex = re.compile(reString)
    
    if url is None:
        url = ('https://rbsp-ect.lanl.gov/data_pub/rbsp{}/MagEphem/'
            'definitive/{}/'.format(sc_id.lower(), splitDate[0]))
        
    response = urllib.request.urlopen(url)
    data = response.read().decode()
    # Use regular expressions to find the filename of intrest.
    mo = magEphemRegex.findall(data)
    assert len(mo) == 1, 'Error, none or multiple files found! {}'.format(mo)
    dataUrl = os.path.join(url, mo[0].split('\"')[0])
    return dataUrl

def saveRBSPMagEphem(sc_id, date, fPath, fName = None, **kwargs):

    url = kwargs.get('url', None)
    Bmodel = kwargs.get('Bmodel', 'T89D')
    fType = kwargs.get('fType', 'txt')
    
    # Create output file directory if necessary
    if not os.path.exists(fPath):
        os.makedirs(fPath)
        print('Created directory:', fPath)
    
    # load and save data
    try:
        dataUrl = findRBSPMagEphemUrl(sc_id, date, **kwargs)
    except AssertionError as err:
        print(str(err))
        return
        
    # If file name not supplied, keep the old one.
    if fName is None:
        fName = dataUrl.split('/')[-1]
    print(fPath, fName)
    urllib.request.urlretrieve(dataUrl, os.path.join(fPath, fName))
    return

if __name__ == '__main__':
    import argparse
    import itertools
    parser = argparse.ArgumentParser(description=('This script downloads the '
            'RBSP magEpehem data from a website defined in the code.'))
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
    parser.add_argument('-m', '--model', default='T89D', help=('Magnetic field '
                                    'model to download the MagEphem files for.'))
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

    delta = endDate - startDate
    for (sc_id, i) in itertools.product(args.sc_id, range(delta.days + 1)):
        d = startDate + timedelta(days=i)
        saveRBSPMagEphem(sc_id, d, '/home/mike/research/rbsp/magephem/'
                'rbsp{}'.format(sc_id.lower()), Bmodel=args.model)


    # date = datetime(2015,2,2)
    # sc = 'b'

    # startDate = datetime(2017, 7, 1)
    # endDate =  datetime(2017, 8, 1)
    # delta = endDate - startDate
    
    # for sc in ['A', 'B']:
    #     for i in range(delta.days + 1):
    #         date = startDate + timedelta(days=i)
    #         print(date)
    #         saveRBSPMagEphem(sc, date, '/home/mike/research/rbsp/magephem/'
    #             'rbsp{}'.format(sc.lower()))
