#!/bin/bash
# Read how to create more here: http://linuxcommand.org/wss0010.php

echo "Downloading FIREBIRD TLEs"
#scp -r mshumko@europa.ssel.montana.edu:/srv/sftp/firebird-project/firebird-data/MSU_Side/TLEs/msu_tles* /home/mike/research/firebird/MSU_Side/tle

rsync -t mshumko@europa.ssel.montana.edu:/srv/sftp/firebird-project/firebird-data/MSU_Side/TLEs/msu_tles* /home/mike/research/firebird/MSU_Side/tle
