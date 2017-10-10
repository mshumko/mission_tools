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
    assert len(mo) == 1, 'Error, none or multiple files found!'
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
    #date = datetime(2017,3,10)
    sc_id = ['A','B']

    startDate = datetime(2017, 7, 1)
    endDate =  datetime(2017, 8, 1)
    delta = endDate - startDate
    
    for sc in sc_id:
        for i in range(delta.days + 1):
            date = startDate + timedelta(days=i)
            print(date)
            saveRBSPMagEphem(sc, date, '/home/mike/research/rbsp/magephem/'
                'rbsp{}'.format(sc.lower()))
