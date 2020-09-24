#download_fb_data_web.sh
#!/bin/bash

echo "Downloading processed FB data from the web."
wget -N -r --accept "FU*" --no-directories --no-parent -e robots=off -P /home/mike/research/firebird/Datafiles/FU_3/context/level2 http://solar.physics.montana.edu/FIREBIRD_II/Data/FU_3/context/
wget -N -r --accept "FU*" --no-directories --no-parent -e robots=off -P /home/mike/research/firebird/Datafiles/FU_4/context/level2 http://solar.physics.montana.edu/FIREBIRD_II/Data/FU_4/context/

wget -N -r --accept "FU*" --no-directories --no-parent -e robots=off -P /home/mike/research/firebird/Datafiles/FU_3/hires/level2 http://solar.physics.montana.edu/FIREBIRD_II/Data/FU_3/hires/
wget -N -r --accept "FU*" --no-directories --no-parent -e robots=off -P /home/mike/research/firebird/Datafiles/FU_4/hires/level2 http://solar.physics.montana.edu/FIREBIRD_II/Data/FU_4/hires/
