# model.py

# HYSPLITm - HYSPLIT Manager

# The HYSPLIT model is maintained by the NOAA Air Resources Lab. Users of this
# program should properly credit and reference NOAA ARL.
# For more information please visit:
#   http://www.arl.noaa.gov/HYSPLIT_info.php
#   http://www.arl.noaa.gov/disclaimer.php

# See the README for more information.


import sys
import os
import pickle
import shutil
from time import sleep
import datetime as dt
import subprocess
import re
import multiprocessing
import numpy as np


class HySetup(object):
    """Sets up the HYSPLIT model.
    """

    def __init__(self, parent, n_exe=-1):
        """
        If n_exe is -1, will query user for number of instances of HYSPLIT
        to run simultaneously.
        """
        self._parent = parent
        # Setup directories
        self.h_dir = self._parent._h_dir
        self.var_dir = self._parent._var_dir
        self.work_dir = self._parent._work_dir

        # Look for previous pickled setup file
        if os.path.isfile(os.path.join(self.var_dir, 'SETUP')):
            s_file = open(os.path.join(self.var_dir, 'SETUP'), 'r')
            self.met_dir, self.out_dir, self.n_exe, \
            self.auto_dl, self.incl_ld = pickle.load(s_file)
            s_file.close()
            # Overwrite number of executables to run if necessary
            if self.n_exe != n_exe and n_exe > 0:
                self.n_exe = n_exe
        else:
            # Get HYSPLIT I/O directories
            self.met_dir = raw_input(
                'Enter HYSPLIT meteorology files directory: ')
            self.out_dir = raw_input(
                'Enter HYSPLIT output files directory: ')

            # Get number of HYSPLIT instances to run
            if n_exe == -1:
                ncpu = multiprocessing.cpu_count()
                print '''
        There are %s processor cores on this computer.
        Enter the number of simultaneous HYSPLIT executables you want.
                ''' % ncpu
                while True:
                    try:
                        run_mode = input()
                        if type(run_mode).__name__ == 'int' \
                                and 0 < run_mode <= ncpu:
                            break
                        else:
                            print 'Invalid Entry'
                    except (ValueError, SyntaxError, NameError):
                        print 'Invalid Entry'
                if run_mode:
                    self.n_exe = run_mode
                else:
                    self.n_exe = multiprocessing.cpu_count() - 1
            elif n_exe > 0:
                self.n_exe = n_exe
            else:
                self.n_exe = multiprocessing.cpu_count() - 1

            # Get automatic met file download preferenceprint
            # NOTE: haven't fully tested, not included in this version
            if False:
                print '''
        Automatically download missing met files from NOAA ARL server?
        1 for yes, 0 for no:
                '''
                while True:
                    try:
                        auto_dl = input()
                        if type(auto_dl).__name__ == 'int' \
                                and 0 <= auto_dl <= 1:
                            break
                        else:
                            print 'Invalid Entry'
                    except (ValueError, SyntaxError, NameError):
                        print 'Invalid Entry'
                self.auto_dl = bool(auto_dl)
            else:
                self.auto_dl = False

            # Ask if a receptor location directory should be included in output
            print '''
        Save HYSPLIT output files in directories for each receptor?
        1 for yes, 0 for no:
            '''
            while True:
                try:
                    incl_ld = input()
                    if type(incl_ld).__name__ == 'int' \
                            and incl_ld >= 0 and incl_ld <= 1:
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            self.incl_ld = bool(incl_ld)

            # Pickle HYSPLIT Setup file
            s_file = open(os.path.join(self.var_dir, 'SETUP'), 'w')
            pickle.dump([self.met_dir, self.out_dir,
                         self.n_exe, self.auto_dl, self.incl_ld], s_file)
            s_file.close()

        # Ensure there are enough HYSPLIT working directories for each instance
        for i in range(self.n_exe):
            if not (os.path.isdir(os.path.join(
                    self.work_dir, 'h%s' % (i + 1)))):
                shutil.copytree(os.path.join(self.work_dir, 'hdefault'),
                                os.path.join(self.work_dir, 'h%s' % (i + 1)))


class HyBase(object):
    """Inheritable class with methods for using HYSPLIT model and its output.
    """

    def __init__(self, parent):
        self._parent = parent

    def check_date(self, date):
        """Checks if HYSPLIT dates requested are valid dates.
        Date input format (as string): 'YYYYMMDD'
        Returns True for valid, False for invalid date
        """
        try:
            dt.datetime.strptime(date, '%Y%m%d')
            return True
        except ValueError:
            return False

    def run_dir(self, run):
        """Returns proper directory for run with trailing slash.
        """
        if not hasattr(self, 'HyS'):
            self.HyS = self._parent.HyS

        if self.HyS.incl_ld:
            run_dir = os.path.join(self.HyS.out_dir, self.get_var(run)['met'],
                                   run[0:14], '')
        else:
            run_dir = os.path.join(self.HyS.out_dir, self.get_var(run)['met'])
            run_dir = os.path.join(run_dir, '')
        return run_dir

    def get_name(self, lat, lon, met, ver, dur, alt_type, alt, yr, mo, day, hr, non_d=None):
        """Returns the proper name for a HYSPLIT output file.
        If improper inputs, returns False.
        non_d will hold non-default run options if called.
        NOTE: the duration variable, dur, has been updated to four digits. This
        version is therefore not reverse compatible with runs from older
        versions.
        Parameters:
            lat:        Latitude; int or float. Negative for South.
            lon:        Longitude; int or float. Negative for West.
            met:        Met Dataset; string.
                options: 'gdas1', 'gdas0p5', 'edas40', 'rean', 'cmp4'
            ver:        Vertical Motion; int. 0 to 5, 8 now?. See HYSPLIT user's guide.
                        0 is typical and will use the model vertical velocity.
                options: 0:data 1:isob 2:isen 3:dens 4:sigma 5:diverg 6:msl2agl 7:average 8:damped
            dur:        Trajectory Duration (hours); int. Negative for backtrajectory.
            alt_type:   Altitude Type; int.
                options: 0:AGL, 1:MSL
            alt:        Altitude (m): int.
            yr:         Year; int.
            mo:         Month; int.
            day:        Day; int.
            hr:         Hour; int.
            non_d:      Dictionary with non-default HYSPLIT run options. Leave
                        as None to use only default options.
        """
        n = ''
        n = n + ('%05.2f' % abs(lat))[0:2] + '-' + ('%05.2f' % abs(lat))[3:5]
        if lat > 0:
            n = n + 'N_'
        else:
            n = n + 'S_'
        n = n + ('%06.2f' % abs(lon))[0:3] + '-' + ('%06.2f' % abs(lon))[4:6]
        if lon > 0:
            n = n + 'E_'
        else:
            n = n + 'W_'
        n = n + met + '_'
        n = n + str(ver)
        if dur < 0:
            n = n + 'B'
        if dur > 0:
            n = n + 'F'
        n = n + '%04d' % abs(dur)
        if alt_type == 0:
            n = n + 'A'
        if alt_type == 1:
            n = n + 'M'
        n = n + '%05d' % alt
        n = n + '%04d' % yr + '%02d' % mo + '%02d' % day + '%02d' % hr
        # FIXME: should just get key value pairs from non_d and append.
        if non_d is not None:
            non_d_list = ['mdom', 'met', 'mgmin', 'tratio', 'delt', 'tout', 'khmax']
            for i in non_d_list:
                try:
                    n = n + '_' + i + str(non_d[i])
                except (KeyError, TypeError):
                    pass
        return n

    def get_var(self, run):
        """Returns a dictionary of variables with Run instructions.
        Approximate opposite of get_name().
        Returns False on bad Run name.
        Variable dictionary for Run:
            'lat','lon','met','ver','dur','alt_type',
            'alt','yr','mo','day','hr','non_d'
        If non-default options are set, will return a dictionary of
        non-default options.
        in the value for 'non_d'. Otherwise 'non_d' has a value of 'None'.
        NOTE: the duration variable, dur, has been updated to four digits. This
        version is therefore not reverse compatible with runs from older
        versions.
        """
        try:
            run = run.replace('-', '.')
            # "R":   0=Lat, 1=Lon, 2=Met, 3=main run options,
            #       (4+)=non-default options if present
            R = run.rsplit('_')
            num_non_d = len(R) - 4
            lat = float(R[0][0:5])
            if R[0][5] == 'S':
                lat = lat * -1.
            lon = float(R[1][0:6])
            if R[1][6] == 'W':
                lon = lon * -1.
            met = R[2]
            ver = int(R[3][0])
            dur = int(R[3][2:6])
            if R[3][1] == 'B':
                dur = dur * -1
            if R[3][6] == 'A':
                alt_type = 0
            else:
                alt_type = 1
            alt = int(R[3][7:12])
            yr = int(R[3][12:16])
            mo = int(R[3][16:18])
            day = int(R[3][18:20])
            hr = int(R[3][20:22])
            # Non-Default run options
            # FIXME: Not sure what the met_lets is for. Remnant from previous
            # naming scheme, but shouldn't be needed any longer.
            # Actually, this whole section is weird, need to redo and simplify.
            non_d = {}
            re_var_name = re.compile('[a-z]*')
            for i in range(num_non_d):
                tmp = re_var_name.match(R[4 + i])
                var_name = tmp.group()
                val = R[4 + i][tmp.end():]
                if var_name == 'mdom':
                    non_d[var_name] = val
                elif var_name == 'met':
                    vals = {}
                    met_vars = ['tm_pres', 'tm_tpot', 'tm_tamb', 'tm_rain',
                                'tm_mixd', 'tm_relh', 'tm_sphu', 'tm_mixr',
                                'tm_dswf', 'tm_terr']
                    met_lets = ['P', 'O', 'A', 'R', 'D', 'H', 'E', 'W', 'S', 'T']
                    for j in range(10):
                        if val.find(met_lets[j]) == -1:
                            vals[met_vars[j]] = 0
                        else:
                            vals[met_vars[j]] = 1
                    non_d[var_name] = vals
                else:
                    non_d[var_name] = val
            if not non_d:
                # FIXME: Not sure why I put z in here - look into later
                z = None  # analysis:ignore
            return {'lat': lat, 'lon': lon, 'met': met, 'ver': ver, 'dur': dur,
                    'alt_type': alt_type, 'alt': alt, 'yr': yr, 'mo': mo,
                    'day': day, 'hr': hr, 'non_d': non_d}
        except (ValueError, IndexError):
            return False

    def check_run(self, run):
        """Checks if HYSPLIT 'run' is a valid name.
        Boolean return.
        """
        if self.get_var(run) is False:
            return False
        return True

    def check_ex(self, run, mkdir=False):
        """Checks if HYSPLIT output file and directory exists.
        If mkdir is True, will create the directory if needed.
        Returns None if bad run name.
        Returns False if valid name and file doesn't exist.
        Returns True if valid name and file does exist.
        """
        if not self.check_run(run):
            return None
        run_dir = self.run_dir(run)
        if not (os.path.isdir(
                os.path.join(self.HyS.out_dir, self.get_var(run)['met']))) \
                and mkdir:
            os.mkdir(os.path.join(self.HyS.out_dir, self.get_var(run)['met']))
        if not (os.path.isdir(run_dir)):
            os.mkdir(run_dir)
        if os.path.isfile(os.path.join(run_dir, run)):
            return True
        else:
            return False

    def _list_runs(self):
        """Original method for adding runs to a list by user input.
        Returns a list of run names, or False if no valid runs to return (need
        to rerun input).
        NOTE: This is a very old method from an earlier version intended to
        add runs for a specific situation. There are much better ways to do this
        now, but it is left in for the time being as a simple way to add runs
        in interactive mode.
        """
        # -Get needed runs-
        print 'Start HYSPLIT Runs for:'
        print ''
        print 'Number of Locations:'
        while True:
            try:
                num_loc = input()
                if isinstance(num_loc, int) and num_loc >= -1:
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        if num_loc == -1:
            return False
        print ''
        print 'Location (Lat, Long to two decimal places as [X.XX,X.XX]):'
        print '  Note that south lat and west long should be negative'
        lat = []
        lon = []
        for i in range(num_loc):
            print 'Location %s:' % (i + 1)
            while True:
                try:
                    lat_lon = input()
                    if isinstance(lat_lon, list) and len(lat_lon) == 2:
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            lat.append(lat_lon[0])
            lon.append(lat_lon[1])
        print ''
        print '0: Gdas1'
        print '1: Gdas0p5'
        print '2: Edas40'
        print '3: Reanalysis'
        print '4: Coamps4'
        print 'Met Dataset:'
        while True:
            try:
                met_int = input()
                if isinstance(met_int, int) and 0 <= met_int < 5:
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        met_dict = {0: 'gdas1', 1: 'gdas0p5', 2: 'edas40', 3: 'rean', 4: 'cmp4'}
        if met_int == 4:
            print 'WARNING: Naming of HYSPLIT COAMPS met datasets are not'
            print 'unique as the grid varies based on the COAMPS run.'
            print 'Care should be taken to ensure the met datasets for the'
            print 'proper run are included in the HYSPLIT met files directory.'
        met = met_dict[met_int]
        print ''
        print 'Vertical Motion option:'
        while True:
            try:
                ver = input()
                if isinstance(ver, int) and 0 <= ver <= 5:
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        print ''
        print 'Duration hours (negative for backtrajectories):'
        while True:
            try:
                dur = input()
                if isinstance(dur, int):
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        print ''
        print 'Altitude Type, 0:AGL 1:MSL:'
        while True:
            try:
                alt_type = input()
                if isinstance(alt_type, int) and 0 <= alt_type <= 1:
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        print ''
        print 'Number of Altitudes:'
        while True:
            try:
                num_alt = input()
                if isinstance(num_alt, int):
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        print ''
        print 'Altitude (meters):'
        alt = []
        for i in range(num_alt):
            print 'Altitude ' + str(i + 1) + ':'
            while True:
                try:
                    tmp = input()
                    if isinstance(tmp, int):
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            alt.append(tmp)
        print ''
        print 'Date Range (Start Date, End Date as [YYYY,MM,DD,YYYY,MM,DD]):'
        while True:
            try:
                date_range = input()
                if isinstance(date_range, list) and len(date_range) == 6:
                    if self.check_date(str(date_range[0]) +
                                               str(date_range[1]) +
                                               str(date_range[2])) == True \
                            and self.check_date(str(date_range[3]) +
                                                        str(date_range[4]) +
                                                        str(date_range[5])) == True:
                        break
                    else:
                        print 'Invalid Date'
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        b_yr, b_mo, b_day, e_yr, e_mo, e_day = date_range
        print ''
        print 'Run at Hours each day: ([X,X,X] etc.) or'
        print '-1 for [0,3,6,9,12,15,18,21]:'
        while True:
            try:
                hrs = input()
                if hrs == -1:
                    hrs = [0, 3, 6, 9, 12, 15, 18, 21]
                    break
                else:
                    if isinstance(hrs, list) and len(hrs) <= 24:
                        break
                    else:
                        print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        print ''
        print 'Change any Default options? (0:No 1:Yes):'
        while True:
            try:
                nd = input()
                if type(nd).__name__ == 'int' and 0 <= nd <= 1:
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        if nd:
            print 'Enter value for Non-Default run options.'
            print 'Enter "-1" to leave as default.'
            non_d = {}
            while True:
                try:
                    met_ = input('met (\'string\'): ')
                    if isinstance(met, str):
                        break
                    else:
                        if type(met).__name__ == 'int' and met == -1:
                            break
                        else:
                            print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            if met_ != -1:
                non_d['met'] = met_
            while True:
                try:
                    mgmin = input('subgrid size (integer): ')
                    if isinstance(mgmin, int) and (mgmin >= 0 or mgmin == -1):
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            if mgmin != -1:
                non_d['mgmin'] = mgmin
            while True:
                try:
                    tratio = input('Grid cell ratio fraction (float): ')
                    if isinstance(tratio, float) and 0. <= tratio <= 1.:
                        break
                    else:
                        if isinstance(tratio, int) and tratio == -1:
                            break
                        else:
                            print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            if tratio != -1:
                non_d['tratio'] = tratio
            if tratio == -1:
                while True:
                    try:
                        delt = input('integration time step (integer): ')
                        if isinstance(delt, int) \
                                and ((0 <= delt <= 60) or delt == -1):
                            break
                        else:
                            print 'Invalid Entry'
                    except (ValueError, SyntaxError, NameError):
                        print 'Invalid Entry'
                if delt != -1:
                    non_d['delt'] = delt
            while True:
                try:
                    tout = input(
                        'Endpoint output frequency in minutes (integer): ')
                    if isinstance(tout, int) and (tout >= 0 or tout == -1):
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            if tout != -1:
                non_d['tout'] = tout
            while True:
                try:
                    khmax = input(
                        'maximum trajectory duration in hours (integer): ')
                    if isinstance(khmax, int) and (khmax >= 0 or khmax == -1):
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            if khmax != -1:
                non_d['khmax'] = khmax
                non_d['mhrs'] = khmax
            while True:
                try:
                    mdom = \
                        input('top of model domain in meters AGL (integer): ')
                    if isinstance(mdom, int) and (mdom >= 0 or mdom == -1):
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            if mdom != -1:
                non_d['mdom'] = mdom
        else:
            non_d = None
        # -Print options and confirm-
        print ''
        print ''
        print 'Current HYSPLIT Run Request'
        print 'Location(s) (Lat(s), Lon(s)): ', lat, lon
        print 'Met Dataset: ', met
        print 'Vertical Motion Option: ', ver
        print 'Trajectory Duration (negative for backtrajectory): ', dur
        if alt_type == 0:
            tmp = ' AGL'
        else:
            tmp = ' MSL'
        print 'Altitude(s): ', alt, tmp
        print 'Date Range: ', date_range
        print 'Hours to Run at: ', hrs
        if nd:
            print 'Non-Default options: ', non_d
        print 'Is this correct? (0:No 1:Yes): '
        while True:
            try:
                correct = input()
                if isinstance(correct, int) and 0 <= correct <= 1:
                    break
                else:
                    print 'Invalid Entry'
            except (ValueError, SyntaxError, NameError):
                print 'Invalid Entry'
        if correct == 0:
            return False
        # -Generate each run time-
        # For n Runs (at one time each) of HYSPLIT there is an index which has
        # each time value at Yr[n], Mo[n] etc.
        n = 0  # The number of times HYSPLIT is to be run for each case
        yr = []
        mo = []
        day = []
        hr = []
        now = dt.date(b_yr, b_mo, b_day)
        end = dt.date(e_yr, e_mo, e_day)
        while now <= end:
            for i in hrs:
                yr.append(now.year)
                mo.append(now.month)
                day.append(now.day)
                hr.append(i)
                n += 1
            now = now + dt.timedelta(1)
        # -Save to list and return-
        ret = []
        for i_loc in range(num_loc):
            for i_alt in range(num_alt):
                for i in range(n):
                    ret.append(self.get_name(lat[i_loc], lon[i_loc], met,
                                             ver, dur, alt_type, alt[i_alt],
                                             yr[i], mo[i], day[i], hr[i],
                                             non_d))
        return ret

    def add_queue(self, ret_runset=False):
        """Queries user to add more runs to HYSPLIT queue.
        Parameters:
            ret_runset:     If True, will return a list of all newly added
                            run names.
                            If False, will return None
        """
        runset = []
        while True:
            # User add runs
            rslt = self._list_runs()
            more = 0
            if rslt:
                runset += rslt
            else:
                print 'Warning: No runs added'
            # Check if more runs are needed
            while True:
                try:
                    more = input('Add more runs? (0:No 1:Yes): ')
                    if isinstance(more, int) and 0 <= more <= 1:
                        break
                    else:
                        print 'Invalid Entry'
                except (ValueError, SyntaxError, NameError):
                    print 'Invalid Entry'
            # Break user add runs loop if done
            if more == 0:
                break
        # Open Queue.txt to append new runs
        q_file = open(os.path.join(self.HyS.var_dir, 'Queue.txt'), 'a')
        for r in runset:
            q_file.write(r + '\n')
        print '%s HYSPLIT runs added to Queue' % len(runset)
        q_file.close()
        # Return runset if requested
        if ret_runset:
            return runset
        else:
            return None

    def imp_hy_out(self, d, f, ret_meta=False):
        """Import a HYSPLIT output file in directory 'd' and file name 'f'.
        Optional Parameters:
            ret_meta:   If True, will return the metadata of start information
        Return:
            data:       Dictionary of hysplit output data
            metadata:   Optional to return a dictionary of hysplit trajectory
                        start info
        """
        lines = open(os.path.join(d, f), 'r').readlines()
        num_met = int(lines[0].strip().split()[0])
        line = lines[num_met + 1].strip().split()
        traj_num = int(line[0])
        # traj_dir = line[1]
        # traj_vert = line[2]
        if traj_num != 1:
            raise Exception  # Can be removed later if need be
        line = lines[num_met + 2].strip().split()
        metadata = {}
        metadata['st_yr'] = int(line[0])
        metadata['st_mo'] = int(line[1])
        metadata['st_day'] = int(line[2])
        metadata['st_hr'] = int(line[3])
        metadata['st_lat'] = float(line[4])
        metadata['st_lon'] = float(line[5])
        metadata['st_hgt'] = float(line[6])
        line = lines[num_met + 3].strip().split()
        metvar_num = int(line[0])
        metvar = line[1:]
        if len(metvar) != metvar_num:
            raise Exception
        header_lines = 1 + num_met + 3
        data = {'yr': [], 'mo': [], 'day': [], 'hr': [], 'minute': [],
                'fort_hour': [], 'age': [], 'lat': [], 'lon': [], 'hgt': []}
        for x in metvar:
            data[x] = []

        for i in range(header_lines, len(lines)):
            line = lines[i].strip().split()
            data['yr'].append(int(line[2]))
            data['mo'].append(int(line[3]))
            data['day'].append(int(line[4]))
            data['hr'].append(int(line[5]))
            data['minute'].append(int(line[6]))
            data['fort_hour'].append(int(line[7]))
            data['age'].append(float(line[8]))
            data['lat'].append(float(line[9]))
            data['lon'].append(float(line[10]))
            data['hgt'].append(float(line[11]))
            for x in range(metvar_num):
                data[metvar[x]].append(float(line[12 + x]))

        datapoints = len(data['yr'])
        for key in data.keys():
            if len(data[key]) != datapoints:
                raise Exception

        if ret_meta:
            return data, metadata
        return data

    def get_hy_runs(self, runs, hy_outdir,
                    metvars=['PRESSURE', 'THETA', 'AIR_TEMP', 'RAINFALL',
                             'MIXDEPTH', 'RELHUMID', 'SPCHUMID', 'H2OMIXRA',
                             'TERR_MSL', 'SUN_FLUX']):
        """Returns an array of endpoint data for each trajectory run in 'runs'.
        Parameters:
            runs:       A list of trajectory names that will be imported.
            hy_outdir:  Directory of HYSPLIT output.
            metvars:    A list of the metvars to import. All by default.
        TODO: Currently requires all trajectories to have the same number of
              met vars. Assumed to be all metvars right now.
        """
        num_traj = len(runs)
        variables = ['yr', 'mo', 'day', 'hr', 'minute',
                     'fort_hour', 'age', 'lat', 'lon', 'hgt']
        variables = variables + metvars
        num_vars = len(variables)
        run_vars = self.get_var(runs[0])
        fdir = os.path.join(hy_outdir, run_vars['met'], runs[0][:14])
        first_traj = self.imp_hy_out(fdir, runs[0])
        num_endpoints = len(first_traj['age'])
        data = np.zeros([num_traj, num_endpoints, num_vars], dtype=float)
        data[:] = np.NaN
        # for each trajectory
        for t in range(num_traj):
            fdir = os.path.join(hy_outdir, run_vars['met'], runs[t][:14])
            traj = self.imp_hy_out(fdir, runs[t])
            # for each endpoint
            for e in range(num_endpoints):
                # for each variable
                for v in range(num_vars):
                    try:
                        data[t, e, v] = traj[variables[v]][e]
                    except IndexError:
                        data[t, e, v] = np.NaN

        return data, variables

    def _dl_met(self, met_names, met_type, met_dir):
        """Downloads requested met files from the NOAA ARL server.
        Returns True if downloads were successful, False otherwise.
        #TODO: not tested yet. Need to add other met types
        Download files manually at:
        ftp://arlftp.arlhq.noaa.gov/archives
        """
        import ftplib
        if isinstance(met_names, str):
            met_names = [met_names]
        try:
            arl = ftplib.FTP(r'arlftp.arlhq.noaa.gov')
            arl.login('anonymous')
            arl.cwd('archives/')
            if met_type is 'gdas1':
                arl.cwd('gdas1')
                for mn in met_names:
                    f = open(os.path.join(met_dir, 'Gdas1', mn), 'wb')
                    arl.retrbinary('RETR %s' % mn, f.write)
                    f.close()
            else:
                raise Exception('met_type not found')
            return True
        except:
            print 'ARL MET DL WARNING for %s:' % mn, sys.exc_info()
            return False


class HyRun(multiprocessing.Process, HyBase):
    """Runs HYSPLIT executable for any requests in Queue.
    """

    def __init__(self, HyS, self_dir, HyQ, ELock, TLock):
        """Sets needed class variables.
        """
        # Class variables
        self.HyS = HyS
        self.self_dir = self_dir  # Working Directory for this HyRun instance
        self.HyQ = HyQ
        self.ELock = ELock
        self.TLock = TLock
        self.error_file = os.path.join(self.HyS.var_dir, 'Error.txt')
        self.trunc_file = os.path.join(self.HyS.var_dir, 'Truncated.txt')
        # Initialize process and run
        super(HyRun, self).__init__()
        self.run()

    def setup_file(self, run, supress_out_var=False):
        """Creates or modifies 'SETUP.CFG' file for run.
        """
        var = self.get_var(run)
        defaults = {
            'tratio': 0.75,
            'delt': 0.0,
            'mgmin': 10,
            'khmax': 9999,
            'kmixd': 0,
            'kmsl': 0,
            'nstr': 0,
            'mhrs': 9999,
            'nver': 0,
            'tout': 60,
            'tm_pres': 1,
            'tm_tpot': 1,
            'tm_tamb': 1,
            'tm_rain': 1,
            'tm_mixd': 1,
            'tm_relh': 1,
            'tm_sphu': 1,
            'tm_mixr': 1,
            'tm_dswf': 1,
            'tm_terr': 1,
            'dxf': 1.00,
            'dyf': 1.00,
            'dzf': 0.01
        }
        if supress_out_var:
            for dn in ['tm_pres', 'tm_tpot', 'tm_tamb', 'tm_rain', 'tm_mixd',
                       'tm_relh', 'tm_sphu', 'tm_mixr', 'tm_dswf', 'tm_terr']:
                defaults[dn] = 0
        w_lst = ['tratio', 'delt', 'mgmin', 'khmax', 'kmixd', 'kmsl', 'nstr', 'mhrs',
                 'nver', 'tout',
                 'tm_pres', 'tm_tpot', 'tm_tamb', 'tm_rain', 'tm_mixd', 'tm_relh',
                 'tm_sphu', 'tm_mixr', 'tm_dswf', 'tm_terr',
                 'dxf', 'dyf', 'dzf']
        s_file = open(os.path.join(self.self_dir, 'SETUP.CFG'), 'w')
        s_file.write(' &SETUP\n')
        for i in w_lst:
            # Because I actually stored the kmsl variable in 'alt_type' rather
            # than 'non_d', look to see if this is the var in w_lst first
            # FIXME: The whole ways I use the 'non_d' dictionary may need
            #   reworking. It's a bit convoluted right now.
            if i is 'kmsl':
                tmp = var['alt_type']
                s_file.write(' %s = %i,\n' % (i, tmp))
            else:
                try:
                    tmp = var['non_d'][i]
                    s_file.write(' ' + i + ' = ' + str(tmp) + ',\n')
                except (KeyError, TypeError):
                    try:
                        tmp = var['non_d']['met'][i]
                        s_file.write(' ' + i + ' = ' + str(tmp) + ',\n')
                    except (KeyError, TypeError):
                        s_file.write(' ' + i + ' = ' + str(defaults[i]) + ',\n')
        s_file.write(' /\n')
        s_file.close()
        return True

    def control_file(self, run):
        """Creates or modifies 'CONTROL' file for run.
        """
        var = self.get_var(run)
        num_met, met_dirs, met_name = self.met_files(run)
        c_file = open(os.path.join(self.self_dir, 'CONTROL'), 'w')
        c_file.write(('%04d' % var['yr'])[2:4] +
                     ' %02d %02d %02d\n' % (var['mo'], var['day'], var['hr'])
                     )
        c_file.write('1\n')
        c_file.write('%s %s %s\n' % (var['lat'], var['lon'], var['alt']))
        c_file.write('%s\n' % var['dur'])
        c_file.write('%s\n' % var['ver'])
        try:
            spam = var['non_d']['mdom']
        except (KeyError, TypeError):
            spam = 10000
        c_file.write('%07.1f\n' % float(spam))
        c_file.write('%s\n' % num_met)
        for i in range(num_met):
            c_file.write(met_dirs[i] + '\n')
            c_file.write(met_name[i] + '\n')
        run_dir = self.run_dir(run)
        c_file.write(run_dir + '\n')
        c_file.write(run + '\n')
        c_file.close()
        return True

    def check_met(self, met_dir, met_name):
        """Checks if requested HYSPLIT met dataset exists.
        If self.HyS.auto_dl is True, will attempt to download any missing
        met files.
        Returns boolean for existance of dataset.
        #TODO: need to add a way to hold execution of these processes while
        downloads complete, but prevent multiple downloads from occuring.
        Probably need a separate process and queue of downloads to add to.
        """
        if os.path.isfile(os.path.join(met_dir, met_name)):
            return True
        else:
            # TODO: not working yet
            return False
            if self.HyS.auto_dl:
                met_type = None  # TODO: get met type (prob pass into method)
                rslt = self._dl_met(met_name, met_type, met_dir)
                if rslt:
                    return True
        return False

    def met_files(self, run):
        """Determines which met files are needed and returns them.
        Returns a list of strings with:
            [number of met files, [list of met dirs],[list of met names]]
        """
        var = self.get_var(run)
        yr = var['yr']
        mo = var['mo']
        sep = os.path.sep
        if var['met'] == 'gdas1':
            # Grabs only met files needed by run
            met_dirs = []
            met_name = []
            Gdas1_dict = {0: 'dec', 1: 'jan', 2: 'feb', 3: 'mar', 4: 'apr', 5: 'may',
                          6: 'jun', 7: 'jul', 8: 'aug', 9: 'sep', 10: 'oct',
                          11: 'nov', 12: 'dec'}
            if var['dur'] < 0:
                end = dt.datetime(var['yr'], var['mo'], var['day'])
                st = end + dt.timedelta(hours=var['dur'])
            else:
                st = dt.datetime(var['yr'], var['mo'], var['day'])
                end = st + dt.timedelta(hours=var['dur'])
            cur = st
            while cur.date() <= end.date():
                if cur.day <= 7:
                    wk = 1
                elif cur.day <= 14:
                    wk = 2
                elif cur.day <= 21:
                    wk = 3
                elif cur.day <= 28:
                    wk = 4
                else:
                    wk = 5
                name = 'gdas1.%s%s.w%1i' % (Gdas1_dict[cur.month], ('%04d' % cur.year)[2:4], wk)
                if name not in met_name:
                    met_name.append(name)
                cur += dt.timedelta(days=1)
            num_met = len(met_name)
            for i in range(num_met):
                met_dirs.append(os.path.join(self.HyS.met_dir, 'Gdas1' + sep))
            if num_met > 12:
                raise Exception('HYSPLIT only accepts a max of 12 met files, need to concat them for this run')
        elif var['met'] == 'gdas0p5':
            # daily files for each day of traj
            met_dirs = []
            met_name = []
            st = dt.datetime(var['yr'], var['mo'], var['day'], var['hr'])
            et = st + dt.timedelta(hours=var['dur'])
            num_met = np.abs((st.date() - et.date()).days) + 1
            for i in range(num_met):
                met_dirs.append(os.path.join(self.HyS.met_dir, 'Gdas0p5' + sep))
            for n in range(num_met):
                mt = st + dt.timedelta(days=n * np.sign(var['dur']))
                met_name.append('%04d%02d%02d_gdas0p5' %
                                (mt.year, mt.month, mt.day))
        elif var['met'] == 'edas40':
            # Always return current month, previous month, next month
            met_dirs = []
            met_name = []
            Edas40_dict = {0: 'dec', 1: 'jan', 2: 'feb', 3: 'mar', 4: 'apr',
                           5: 'may', 6: 'jun', 7: 'jul', 8: 'aug', 9: 'sep',
                           10: 'oct', 11: 'nov', 12: 'dec'}
            fmt_yr = yr - 2000
            num_met = 6
            for i in range(num_met):
                met_dirs.append(os.path.join(self.HyS.met_dir, 'Edas40' + sep))
            if mo == 1:
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[12], fmt_yr - 1, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[12], fmt_yr - 1, 2))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo], fmt_yr, 2))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo + 1], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo + 1], fmt_yr, 2))
            elif mo == 12:
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo - 1], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo - 1], fmt_yr, 2))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo], fmt_yr, 2))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[1], fmt_yr + 1, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[1], fmt_yr + 1, 2))
            else:
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo - 1], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo - 1], fmt_yr, 2))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo], fmt_yr, 2))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo + 1], fmt_yr, 1))
                met_name.append('edas.%s%02i.%03i'
                                % (Edas40_dict[mo + 1], fmt_yr, 2))
        elif var['met'] == 'rean':
            # Always return current month, previous month, next month
            met_dirs = []
            met_name = []
            num_met = 3
            for i in range(num_met):
                met_dirs.append(os.path.join(self.HyS.met_dir, 'Reanalysis' + sep))
            if mo == 1:
                met_name.append('RP%04d%02d.gbl' % (yr - 1, 12))
                met_name.append('RP%04d%02d.gbl' % (yr, mo))
                met_name.append('RP%04d%02d.gbl' % (yr, mo + 1))
            elif mo == 12:
                met_name.append('RP%04d%02d.gbl' % (yr, mo - 1))
                met_name.append('RP%04d%02d.gbl' % (yr, mo))
                met_name.append('RP%04d%02d.gbl' % (yr + 1, 1))
            else:
                met_name.append('RP%04d%02d.gbl' % (yr, mo - 1))
                met_name.append('RP%04d%02d.gbl' % (yr, mo))
                met_name.append('RP%04d%02d.gbl' % (yr, mo + 1))
        elif var['met'] == 'cmp4':
            # daily files for each day of traj
            # nominally two grids for now, can add check for more later
            met_dirs = []
            met_name = []
            num_grid = 2
            st = dt.datetime(var['yr'], var['mo'], var['day'], var['hr'])
            et = st + dt.timedelta(hours=var['dur'])
            num_days = np.abs((st.date() - et.date()).days) + 1
            num_met = num_days * num_grid
            for i in range(num_met):
                met_dirs.append(os.path.join(self.HyS.met_dir, 'cmp4' + sep))
            for n in range(num_days):
                ct = st + dt.timedelta(days=n * np.sign(var['dur']))
                for g in range(num_grid):
                    met_name.append('%04d%02d%02d_cmp4_g%d' %
                                    (ct.year, ct.month, ct.day, g + 1))
        else:
            return False, False, False
        # Return met files
        return [num_met, met_dirs, met_name]

    def check_file(self, run):
        """Checks if a HYSPLIT output file was created properly.
        Checks if lines need to be reformatted to longer column width, and
        corrects if needed.
        Checks for the proper number of lines in a file.
        Parameters:
            True if file is now correct
            False if it doesn't exist or no endpoints
            None if file was truncated by HYSPLIT before reaching end of
                expected trajectory duration
        #TODO: May want to do more thorough checking of file validity
        """
        var = self.get_var(run)
        # check file exists
        if self.check_ex(run) is False:
            return False
        sleep(0.1)  # FIXME: Not sure why this is here... file I/O?
        # list of file lines
        l = open(os.path.join(self.run_dir(run), run), 'r').readlines()
        num = len(l)  # number of lines
        num_met = int(l[0].split()[0])  # number of met files (lines)
        h_num = 4 + num_met  # number of header line assuming 1 traj
        e_num = num - h_num  # number of endpoints lines
        # return False if no endpoints
        if e_num <= 0:
            return False
        le1 = len(l[h_num].split())  # number of entries in 1st endpt line
        le2 = len(l[h_num + 1].split())  # number of entries in 2nd endpt line
        if le1 != le2:  # TODO: fix line wrapping needed?
            for s in range(h_num, h_num + (e_num / 2)):
                l[s] = l[s].rstrip() + l.pop(s + 1)
            # rewrite correct file
            open(os.path.join(self.run_dir(run), run), 'w').writelines(l)
            e_num /= 2
        try:
            dur_mult = 60.0 / float(var['non_d']['tout'])
        except (KeyError, TypeError):
            dur_mult = 1
        if e_num < float(var['dur']) * dur_mult:
            # trajectory truncated
            return None
        # return True if all pass
        return True

    def ex_run(self, run):
        """Executes a run as a subprocess with proper setup and control files.
        Also checks for Coamps4 run, which, for now, has all output variables
        turned off.
        Assumes that if running on a posix system, needs to execute via wine
        """
        var = self.get_var(run)
        if var['met'] == 'cmp4':
            sov = True
        else:
            sov = False
        self.setup_file(run, supress_out_var=sov)
        self.control_file(run)
        osname = os.name
        if osname == 'nt':
            precall = ''
        elif osname == 'posix':
            precall = 'wine32 '
        else:
            raise Exception('HYSPLITm only setup for windows and posix/wine systems')
        subprocess.call([precall + os.path.join(self.self_dir, 'hyts_std.exe')],
                        shell=True)

    def run(self):
        """Runs HYSPLIT executable for any requests in Queue file.
        """
        os.chdir(self.self_dir)
        print 'HyRun %s starting from directory %s\n' % (self.name,
                                                         str(self.self_dir))
        # Get system arguments to check for run override
        #   If the -o flag was given, run even though a file exists
        if '-o' in sys.argv:
            override = True
        else:
            override = False
        # Loop to get Run from HyQ until a 'STOP' is found
        for run in iter(self.HyQ.get, 'STOP'):
            # Check that run has valid file name and if it already exists
            rslt = self.check_ex(run, mkdir=True)
            # Determine if a run is needed
            ex_run = False
            conf = False
            if rslt is None:
                # Bad File name
                with self.ELock:
                    open(self.error_file, 'a').write('%s - BADFILENAME\n' % run)
            elif rslt is False:
                # Run doesn't exist, execute run
                ex_run = True
            elif rslt is True:
                if override:
                    # Override any existing file
                    ex_run = True
                else:
                    # File exists, check if correct format
                    conf = self.check_file(run)
                    if conf is False:
                        ex_run = True
            # If run is needed, check for met files, then execute run
            # Executing run sets up SETUP.CFG and CONTROL files and runs
            # HYSPLIT executable.
            if ex_run:
                # Check for existance of needed met files
                num_met, met_dirs, met_name = self.met_files(run)
                met_ex = True
                for i in range(num_met):
                    if not self.check_met(met_dirs[i], met_name[i]):
                        met_ex = False
                if not met_ex:
                    # No Met File
                    with self.ELock:
                        open(self.error_file, 'a').write('%s - NOMET\n' % run)
                else:
                    self.ex_run(run)
                    conf = self.check_file(run)
                    # if conf is False or None, try one more time
                    if not conf:
                        self.ex_run(run)
                        conf = self.check_file(run)
                    # If conf is now not True, save to error or truc log
                    if not conf:
                        if conf is False:
                            # Error in run
                            with self.ELock:
                                open(self.error_file, 'a').write('%s\n' % run)
                        elif conf is None:
                            # Truncated run
                            with self.TLock:
                                open(self.trunc_file, 'a').write('%s\n' % run)
                        else:
                            raise Exception('You should not be here.')
            self.HyQ.task_done()
        return


class HyFor(object):
    """A subset of HYSPLIT builtin Fortran routines translated to python.
    """

    def eqvlat(self, xlat1, xlat2):
        """eqvlat.f
        from fortran file doc:
        "WRITTEN ON 3/31/94 BY Dr. Albion Taylor  NOAA / OAR / ARL"
        Note: Not sure exactly what this function is supposed to do, but
        it looks like there are probably numpy functions that could
        do it better.
        """
        sinl1 = np.sin(np.deg2rad(xlat1))
        sinl2 = np.sin(np.deg2rad(xlat2))
        if np.abs(sinl1 - sinl2) > 0.001:
            al1 = np.log((1. - sinl1) / (1. - sinl2))
            al2 = np.log((1. + sinl1) / (1. + sinl2))
        else:
            tau = (-(sinl1 - sinl2) / (2 - sinl1 - sinl2)) ** 2
            al1 = (2. / (2. - sinl1 - sinl2) *
                   (1. + tau * (1. / 3. + tau * (1. / 5. + tau * (1. / 7.)))))
            tau = ((sinl1 - sinl2) / (2. + sinl1 + sinl2)) ** 2
            al2 = (-2. / (2. + sinl1 + sinl2) *
                   (1. + tau * (1. / 3. + tau * (1. / 5. + tau * (1. / 7.)))))
        return np.rad2deg(np.arcsin((al1 + al2) / (al1 - al2)))
