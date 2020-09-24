import downloadKp
import argparse
import os 

default_save_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ),
            'indicies'))

parser = argparse.ArgumentParser()
parser.add_argument("year", help = 'year(s) to download kp from', nargs = '+')
parser.add_argument("-save_dir", help = 'Save directory for the downloaded kp', 
    default = default_save_dir)
parser.add_argument
args = parser.parse_args()

print(args)
for i in args.year: 
    fname = '{}_kp.txt'.format(i)
    # Download the kp from a year specified by the user flag.
    downloadKp.writeToJSON(i, args.save_dir, fname, newDataFlag=1)
