#!/bin/bash

echo "Downloading processed RBSPA level 3 emfisis data from the web."
wget -r --accept "*4sec-geo*" --no-parent robots=off -P /home/mike/RBSP/data/emfisis/ http://emfisis.physics.uiowa.edu/Flight/RBSP-A/L3/2015/

wget -r --accept "*4sec-geo*" --no-parent robots=off -P /home/mike/RBSP/data/emfisis/ http://emfisis.physics.uiowa.edu/Flight/RBSP-A/L3/2016/

wget -r --accept "*4sec-geo*" --no-parent robots=off -P /home/mike/RBSP/data/emfisis/ http://emfisis.physics.uiowa.edu/Flight/RBSP-A/L3/2017/

# Now move the folders out of their individual directories
#find /home/mike/RBSP/data/emfisis/emfisis.physics.uiowa.edu -name "*.cdf" -exec cp -t /home/mike/RBSP/data/emfisis/ {} +

#rm -r /home/mike/RBSP/data/emfisis/emfisis.physics.uiowa.edu

echo "Downloading processed RBSPB level 3 emfisis data from the web."
wget -r --accept "*4sec-geo*" --no-parent robots=off -P /home/mike/RBSP/data/emfisis/ http://emfisis.physics.uiowa.edu/Flight/RBSP-B/L3/2015/

wget -r --accept "*4sec-geo*" --no-parent robots=off -P /home/mike/RBSP/data/emfisis/ http://emfisis.physics.uiowa.edu/Flight/RBSP-B/L3/2016/

wget -r --accept "*4sec-geo*" --no-parent robots=off -P /home/mike/RBSP/data/emfisis/ http://emfisis.physics.uiowa.edu/Flight/RBSP-B/L3/2017/

# Now move the folders out of their individual directories
find /home/mike/RBSP/data/emfisis/emfisis.physics.uiowa.edu -name "*.cdf" -exec cp -t /home/mike/RBSP/data/emfisis/ {} +

rm -r /home/mike/RBSP/data/emfisis/emfisis.physics.uiowa.edu
