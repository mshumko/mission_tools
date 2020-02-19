# Supported Missions
This git repo holds the collection of Pyhton and shell scripts that download,
read, and visualize data from a handful of missions including:
- AeroCube-6 (AC6)
- The Van Allen Probes (RBSP) - supported instruments are EMFISIS, MagEIS, and RBSPICE
- Geostationary Operational Environmental Satellite (GOES)
- Focused Investigations of Relativistic Electron Burst Intensity, Range, and Dynamics (FIREBIRD-II).

# Dependencies
The major dependencies are listed in requirements.txt file. The matplotlib version
is a werid dependence. If you use the specified version, the map features in the
AC6 and FIREBIRD-II plotting programs will work. Newer versions of matplotlib 
break Cartopy.

# The ```firebird/''' folder contains:
- Download the raw data from the Europa server ```download_fb_data_raw.sh```
- Download the processed data from the web ```download_fb_data_web.sh```
- Download the two line elements (TLEs) ```download_fb_tle.sh```
- Programs that plot the context and hires time series and ground track are ```plot_context.py``` and ```plot_hires.py``` and both have a command line interface
- A ```flag_dropouts.py``` program that can flag times when the HiRes data contains dropouts.

# The ```ac6/``` folder contains:
- A ```download_ac6_data.py``` Python program to download AC6 data from the RBSP Science Gateway. It has a command line interface
- A AC6 data reader as well as a time series and map plotter in ```read_ac_data.py```. Also has a command line interface.

# The ```orbit/``` folder that contains scripts to propagate orbits with the SGP4 algorithm. 
- To run for a mission, save the mission's TLEs in some folder
- Then edit ```make_ephem.py``` and include a path to that TLE folder
- Run ```python3 make_ephem.py -h``` to get help on command line arguments. An example is ```python3 make_ephem.py FU3 FU4 2020 1 1 2020 2 15``` to propagate FU3 and FU4 orbits from 1 January 2020 to 15 Febrary 2020 with the one minute default cadence
- Output ephemeris will be saved in the ```data/``` folder

# The ```misc/``` folder contains:
- 
