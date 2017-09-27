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
        url = ('http://emfisis.physics.uiowa.edu/Flight/RBSP-' + sc_id + '/L2/' + 
        splitDate[0] + '/' + splitDate[1] + '/' +  splitDate[2]) 
          
    response = urllib.request.urlopen(url)
    data = response.read().decode()
    print(hour)
    print('{}-{}_emfisis-L2_{}{}'.format(detectorType, spectraType, 
        ''.join(splitDate), hour))
    ind = data.find('{}-{}_emfisis-L2_{}{}'.format(detectorType, spectraType, 
        ''.join(splitDate), hour))
    lastInd = data.find('.cdf', ind)
    dataUrl = data[ind - 68:lastInd + 4]
    return dataUrl
    
def saveEMFISISSpectra(sc_id, date, fPath, fName = None, **kwargs):

    detectorType = kwargs.get('detectorType', 'WFR')
    spectraType = kwargs.get('spectraType', 'spectral-matrix-diagonal')
    hour = kwargs.get('hour', '')
    
    # Create output file directory if necessary
    if not os.path.exists(fPath):
        os.makedirs(fPath)
        print('Created directory:', fPath)
    
    # load and save data
    dataUrl = findEMFISISSpectraUrl(sc_id, date, spectraType = spectraType, 
    detectorType = detectorType, hour = hour)
    #print(dataUrl.split('/'))
    # If file name not supplied, keep the old one.
    if fName is None:
        fName = dataUrl.split('/')[-1]
    print(fPath, fName)
    urllib.request.urlretrieve(dataUrl, os.path.join(fPath, fName))
    return

if __name__ == '__main__':
    #date = datetime(2017,3,10)
    #sc_id = 'B'
    #dataUrl = findEMFISISSpectraUrl(sc_id, date)
    #saveEMFISISSpectra(sc_id, date, '/home/mike/temp')
    
    startDate = datetime(2014, 6, 1)
    endDate =  datetime(2017, 6, 1)
    delta = endDate - startDate
    for sc_id in ['A', 'B']:
        for i in range(delta.days + 1):
            date = startDate + timedelta(days=i)
            print(date)
            try:
                saveEMFISISSpectra(sc_id, date, '/home/ms30715/ssd_data/rbsp/data'
                '/emfisis/WFR/RBSP{}/L2'.format(sc_id.upper()), 
                    spectraType = 'spectral-matrix-diagonal')
            except urllib.error.HTTPError as err:
                print(str(err))
                continue
