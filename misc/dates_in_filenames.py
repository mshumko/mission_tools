import re, os, glob
import dateutil.parser
import numpy as np

def find_dates(fDir, wildCard = '*', dateTimeConvert = False):
    """
    This will find date instances in a directory with a wildcard
    """
    # List full paths to files in fDir that satisfy the wildcard kwarg
    paths = sorted(glob.glob(os.path.join(fDir, wildCard)))
    dateRegex = re.compile(r'\d{2,4}.\d{1,2}.\d{1,2}')
    
    dates = [dateRegex.findall(i) for i in paths]
    
    # If dateTimeConvert kwarg is True, convert to datetime objects
    if dateTimeConvert:
        for i in range(len(dates)):
            for j in range(len(dates[i])):
                dates[i][j] = dates[i][j].replace('_', '-').replace('.', '-')
                #print(d)
                # Try to convert it.
                try:
                    dates[i][j] = dateutil.parser.parse(dates[i][j])       
                except ValueError:
                    print('Warning string ', dates[i][j], 
                    ' is not in proper format!')
                    dates[i][j] = -9999    
    return paths, np.array(dates)
    
    
if __name__ == '__main__':
    #fDir = '/home/mike/FIREBIRD/Datafiles/FU_3/magephem/'
    fDir = '/home/mike/Dropbox/0_firebird_research/orbit_propagator/TLEs'
    paths, dates = find_dates(fDir, dateTimeConvert = True)

    for i in range(len(paths)):
        print(os.path.basename(paths[i]), dates[i])
