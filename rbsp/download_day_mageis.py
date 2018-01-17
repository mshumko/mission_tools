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
        print('Downloading RBSP-{} MagEIS data from {}'.format(rb_id.upper(), date))
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
            
        # Download data
        urllib.request.urlretrieve(url + fName, os.path.join(fPath, fName))
    return 
    
def findFilename(rb_id, joinDate, data):
    """
    This function uses regular expressions to find the correct
    filename to download
    """
    regex = re.compile(r'rbsp{}_rel\d\d_ect-mageis-L\d_{}_v'
        '\d.\d.\d.cdf'.format(rb_id, joinDate))
    mo = regex.search(data)
    return mo.group()
    
if __name__ == '__main__':
    START_DATE = datetime(2017, 11, 19)
    END_DATE = datetime(2017, 11, 20)
    dates = [START_DATE + timedelta(days=d) for d in 
        range((END_DATE - START_DATE).days)]
    #for (rb_id, day) in itertools.product(['a', 'b'], dates):
    for rb_id in ['a', 'b']:
        fPath = '/home/mike/research/rbsp/data/mageis/rbsp{}'.format(rb_id)
        downloadMagEIS(rb_id, dates, fPath, dataLevel=3)
    

