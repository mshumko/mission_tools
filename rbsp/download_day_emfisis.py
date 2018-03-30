# Download EMFISIS data on date
import urllib.request, os
from datetime import datetime
from datetime import timedelta

#http://emfisis.physics.uiowa.edu/Flight/RBSP-A/L2/2015/02/02/rbsp-a_WFR-waveform-continuous-burst_emfisis-L2_20150202T00_v1.4.3.cdf

def findEMFISISSpectraUrl(sc_id, date, **kwargs):
    """
    This function returns the ftp url of the emfisis Level 2 HFR/WFR spectra.
    spectraType can be spectral-matrix-diagonal, 'spectra', 'spectra-burst'.
    
    If downloading burst data, specify hour as 'T##' where ## is the numerical hour.
    
    sc_id is either 'A' or 'B'
    """
    url = kwargs.get('url', None)
    detectorType = kwargs.get('detectorType', 'WFR')
    spectraType = kwargs.get('spectraType', 'spectral-matrix-diagonal')
    hour = kwargs.get('hour', '') 
    
    splitDate = date.date().isoformat().split('-')
    
    if url is None:
        url = ('http://emfisis.physics.uiowa.edu/Flight/RBSP-' + sc_id.upper() + '/L2/' + 
        splitDate[0] + '/' + splitDate[1] + '/' +  splitDate[2]) 
          
    response = urllib.request.urlopen(url)
    data = response.read().decode()
    ind = data.find('{}-{}_emfisis-L2_{}{}'.format(detectorType, spectraType, 
        ''.join(splitDate), hour))
    lastInd = data.find('.cdf', ind)
    dataUrl = data[ind - 7:lastInd + 4] # For the old website... ind - 68:lastInd + 4
    return url + '/' + dataUrl
    
def saveEMFISISSpectra(sc_id, date, fPath, fName = None, **kwargs):

    detectorType = kwargs.get('detectorType', 'WFR')
    spectraType = kwargs.get('spectraType', 'spectral-matrix-diagonal')
    hour = kwargs.get('hour', '')
    
    # Create output file directory if necessary
    if not os.path.exists(fPath):
        os.makedirs(fPath)
        print('Created directory:', fPath)
    
    # load and save data
    dataUrl = findEMFISISSpectraUrl(sc_id, date, spectraType=spectraType, 
        detectorType=detectorType, hour=hour)
    #print(dataUrl.split('/'))
    # If file name not supplied, keep the old one.
    if fName is None:
        fName = dataUrl.split('/')[-1]
    print('\n', dataUrl, '\n')
    print(fPath, fName)
    urllib.request.urlretrieve(dataUrl, os.path.join(fPath, fName))
    return

if __name__ == '__main__':
    # This section of the code sets up the argument parser.
    import argparse
    import itertools
    parser = argparse.ArgumentParser(description=('This script downloads the '
            'EMFISIS data from a website defined in the code.'))
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

    print('Downloading data from', startDate.date(), 'to', endDate.date())
    delta = endDate - startDate
    for (sc_id, i) in itertools.product(args.sc_id, range(delta.days + 1)):
        d = startDate + timedelta(days=i)
        print('Downloading EMFISIS from: RBSP-{} on {}'.format(sc_id.upper(), d.date()))
        try:
            saveEMFISISSpectra(sc_id, d, '/home/mike/research/rbsp/'
            'data/emfisis/rbsp{}'.format(sc_id.lower()), 
                spectraType = 'spectral-matrix-diagonal')
        except urllib.error.HTTPError as err:
            print(str(err))
            print('RBSP-{}, date={}'.format(sc_id.upper(), d.date()))
            continue


# for sc_id in args.sc_id:
#     for i in range(delta.days + 1):
#         print(sc_id, startDate + timedelta(days=i))

    
    # startDate = datetime(2015, 2, 2)
    # endDate =  datetime(2015, 2, 2)
    # delta = endDate - startDate
    # for sc_id in ['A', 'B']:
    #     for i in range(delta.days + 1):
    #         date = startDate + timedelta(days=i)
    #         print(date)
    #         try:
    #             saveEMFISISSpectra(sc_id, date, '/home/mike/research/rbsp/'
    #             'data/emfisis/rbsp{}'.format(sc_id.lower()), 
    #                 spectraType = 'spectral-matrix-diagonal')
    #         except urllib.error.HTTPError as err:
    #             print(str(err))
    #             continue
