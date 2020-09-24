#!/bin/bash

echo "Downloading raw FIREBIRD data"

echo "FU3"
rsync -tr mshumko@europa.ssel.montana.edu:/srv/sftp/firebird-project/firebird-data/MSU_Side/FU_3/Downlinked_Data/ /home/mike/research/firebird/MSU_Side/FU_3/Downlinked_Data

echo "FU4"
rsync -tr mshumko@europa.ssel.montana.edu:/srv/sftp/firebird-project/firebird-data/MSU_Side/FU_4/Downlinked_Data/ /home/mike/research/firebird/MSU_Side/FU_4/Downlinked_Data

