HYSPLITm - HYSPLIT Manager

The HYSPLIT model is maintained by the NOAA Air Resources Lab. Users of this
program should properly credit and reference NOAA ARL.
For more information please visit:
  http://www.arl.noaa.gov/HYSPLIT_info.php
  http://www.arl.noaa.gov/disclaimer.php

------------------------------
HYSPLIT Manager

Author:   Sam Atwood
          Department of Atmospheric Science
          Colorado State University
email:    satwood@atmos.colostate.edu

This program runs the HYSPLIT trajectory model and manages output files. It
is intended to be used with other analysis scripts which can request specific
HYSPLIT trajectories using a simple, but nominally complete, trajectory
naming system.

This project is licensed under the terms of the MIT license.

Bug reports can be sent to the above email or through github.


USAGE INSTRUCTIONS:
All features of this program are intended to be accessed through the Manager
class. An instance of this class can be created by:
  import HYSPLITm
  H = HYSPLITm.Manager()
HYSPLIT trajectories can be calculated by:
  H.hysplit()
If requested runs are already in the HYSPLIT queue in HYSPLITm/var/Queue.txt,
the runs can be calculated by:
  H.start()
To reset the location of HYSPLIT directories and run setup information, run:
  H.reset()
The last accessed set of runs can be plotted by:
  H.plot()
The meteorological preprocessor for COAMPS4 met files is not included in
this version. Contact Sam Atwood at satwood@atmos.colostate.edu regarding
use of the COAMPS4 met preprocessor.


FIRST RUN INSTRUCTIONS:
This program relies on the hyts_std.exe trajectory model that is part of
the NOAA ARL HYSPLIT program. As this executable is subject to license
agreements with NOAA ARL, it is not distributed with this package. In order
to use this package you must first visit
http://ready.arl.noaa.gov/HYSPLIT.php
First, follow the instructions to download and install the HYSPLIT program. 
Next, find the "hyts_std.exe" executable file from the installed "exec"
directory, e.g. "C:/hysplit4/exec/".
Copy the hyts_std.exe executable file into the "HYSPLITm/working/hdefault/"
directory.
During the first run of HYSPLITm, the program will ask you to configure
the location of several directories the program needs.
    Meteorology Directory: The folder where HYSPLIT met files are stored.
    Output Directory: The folder where HYSPLIT output files will be stored.
Lastly, the preprocessed HYSPLIT met files needed to run hysplit 
(and stored in the Meteorology Directory) can be downloaded from NOAA ARL
via FTP. For more information visit:
http://ready.arl.noaa.gov/archives.php

