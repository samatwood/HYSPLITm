# manager.py

# HYSPLITm - HYSPLIT Manager

# The HYSPLIT model is maintained by the NOAA Air Resources Lab. Users of this
# program should properly credit and reference NOAA ARL.
# For more information please visit:
#   http://www.arl.noaa.gov/HYSPLIT_info.php
#   http://www.arl.noaa.gov/disclaimer.php

# See the README for more information.


# ---Modules---
import sys
import os
import multiprocessing
# import numpy as np
import HYSPLITm.model as model
import HYSPLITm.plot as hp


# import HYSPLITm.met as hm


class Manager(object):
    """The manager that controls all functions of HYSPLITm.
    Class Attributes:
        _runset:    A list of HYSPLIT run names. Used to carry around the last
                    requested runs for use in other methods or classes.
    """

    def __init__(self, supress_msg=False):
        """Initializes an instance of the HYSPLITm Manager Class.
        """
        # HYSPLITm Startup
        if not supress_msg:
            print '''
        ------------------------------
        HYSPLITm - HYSPLIT Manager - v%s

        Check the README for usage details.
        This program provided without warranty.

        All users of this program should properly credit NOAA ARL.
        For more information, visit
        http://www.arl.noaa.gov/HYSPLIT_info.php
        ''' % self.__version__
        self._dir()
        # Create a variable that holds the last run set used
        self._runset = None
        # Create instances of HyBase and Plot for use in interactive mode
        # FIXME: HyB and HyS usage is too complicated. I think due to trying to
        #        both compose and inherit them.
        self.HyB = model.HyBase(self)
        self.HyS = model.HySetup(self)
        self.Plot = hp.Plot(self)
        self.Plot.HyS = model.HySetup(self)

    def _dir(self):
        """Gets directories associated with HYSPLITm.
        """
        self._h_dir = os.path.abspath(os.path.dirname(__file__))
        self._var_dir = os.path.join(self._h_dir, 'var')
        self._work_dir = os.path.join(self._h_dir, 'working')
        self._lib_dir = os.path.join(self._h_dir, 'lib')
        self._hylib_dir = os.path.join(self._lib_dir, 'hylib')
        self._hyutil_dir = os.path.join(self._lib_dir, 'hyutil')
        self._exe_dir = os.path.join(self._h_dir, 'exe')
        if not os.path.isdir(self._var_dir):
            os.mkdir(self._var_dir)
        if not os.path.isdir(self._work_dir):
            raise Exception('No HYSPLIT working directory found')
        for i in ['Error.txt', 'Queue.txt', 'Truncated.txt']:
            if not os.path.isfile(os.path.join(self._var_dir, i)):
                open(os.path.join(self._var_dir, i), 'w').close()

    def reset(self):
        """Resets setup files, including run prefs and directory locations.
        """
        if os.path.isfile(os.path.join(self._var_dir, 'SETUP')):
            os.remove(os.path.join(self._var_dir, 'SETUP'))

    # ---HYSPLIT model methods---
    def start(self, add_runs=False, go=True):
        """A simple start function to call the hysplit() run method.
        This is basically just to allow the HYSPLITm module to be imported and
        run by other scripts or in python interactive mode.
        Other features can be added to this function later if need be.
        Parameters:
            add_runs:   False by default, no runs added to HYSPLIT queue.
            go:         True by default starts running the HYSPLIT model with
                        runs already in the queue.
        """
        self.hysplit(add_runs, go)

    def hysplit(self, add_runs=True, go=True):
        """Main function to run the HYSPLIT model.
        This method uses the /var/Queue.txt file to store run names so runs can
        be added from anywhere, including other programs, and executed by this
        method.
        Parameters:
            add_runs:   If True (Default), will add runs to HYSPLIT queue.
                        Otherwise, will just use HYSPLIT runs already in
                        the queue.
            go:         If True (Default), will start exectuing runs from the
                        queue. If False, will just setup files so the HYSPLIT
                        model is ready, and can add runs.
        """

        ELock = multiprocessing.Lock()
        TLock = multiprocessing.Lock()

        # Setup HYSPLIT Model
        HyS = model.HySetup(self)
        HyB = model.HyBase(self)
        HyB.HyS = HyS

        # Add new runs to HYSPLIT Queue if requested
        if add_runs:
            self._runset = HyB.add_queue(ret_runset=True)

        if go:
            # Start HyQueue to manage Queue
            HyQ = multiprocessing.JoinableQueue(30)
            Qmgr = multiprocessing.Process(target=q_mgr,
                                           args=(HyS, HyQ))
            Qmgr.start()

            # Start the requested number of HyRun instances
            HyRunInst = []
            for i in range(HyS.n_exe):
                print 'Starting Process %s' % (i + 1)
                p = multiprocessing.Process(
                    target=model.HyRun,
                    args=(HyS, os.path.join(HyS.work_dir, 'h%s' % (i + 1)),
                          HyQ, ELock, TLock)
                )
                HyRunInst.append(p)
                p.start()

            # Join Queue Manager until HyQ is empty
            Qmgr.join()
            print 'End of Queue, HYSPLIT Manager shutting down.'

            # Send 'STOP' to HyS.n_exe child processes
            print 'Ending HYSPLIT processes'
            for i in range(HyS.n_exe):
                HyQ.put('STOP')
            for i in range(HyS.n_exe):
                HyRunInst[i].join()
                print 'Joined process %s' % (i + 1)
            print 'Runs Complete'
            print 'HYSPLIT processes shutdown'

        print 'HYSPLITm shutdown'

    # ---Plotting methods---
    def plot(self, runs=None, color_type='date', daily_ticks=False,
             extents=None, lat_space=15, lon_space=25,
             figsize=(14, 8), cax=None, **kwargs):
        """Creates a new HYSPLIT plot of a specific type.
        On first call, will create an instance of the Plot class at self.Plot
        that includes self.Plot.fig and self.Plot.map with the current figure
        and map being plotted.
        Optional Parameters:
            color_type:     The type of plot to run.
                Options:
                    'num':  Color trajectories by trajectory number
                    'date': Color trajectories by date
                    'alt':  Color trajectories by altitude
                    'age':  Color trajectories by age
            extents:        The extents of the map, if needed.
                            Leave as None to autocalculate.
        #TODO:
        Need to include other types of plots here, e.g. alt-lat axes.
        """
        # If a set of runs to plot isn't given, use the last active runset
        if runs is None:
            runs = self._runset
        else:
            self._runset = runs
        # Set up new plot figure
        self.Plot.fig, self.Plot.m = self.Plot.map_setup(extents, figsize,
                                                         r=runs[0], cax=cax,
                                                         lat_space=lat_space, lon_space=lon_space)
        # Plot trajectories
        ret = self.Plot.traj_plot(self.Plot.m, runs,
                                  color_type=color_type, daily_ticks=daily_ticks,
                                  **kwargs)
        if 'cb' in kwargs:
            if kwargs['cb'] == -1:
                return ret

    def rta_plot(self, runs=None, extents=None, lat_space=15, lon_space=25,
                 figsize=(14, 8), cax=None, **kwargs):
        """Creates a new HYSPLIT contoured Residence Time Analysis plot.
        NOTE: This method not fully tested yet, use with caution.
        #TODO:
        """
        # If a set of runs to plot isn't given, use the last active runset
        if runs is None:
            runs = self._runset
        else:
            self._runset = runs
        # Set up new plot figure
        self.Plot.fig, self.Plot.m = self.Plot.map_setup(extents, figsize,
                                                         r=runs[0], cax=cax,
                                                         lat_space=lat_space, lon_space=lon_space)
        # Plot trajectories
        ret = self.Plot.rta_plot(self.Plot.m, runs, **kwargs)

    # ---Met processing methods---
    def met_preprocess(self, met, vert_sys, met_dir):
        """Initializes an instance of the Met() class as self.Met.
        This will allow for conversion of meteorological files into the
        proper ARL format for running HYSPLIT.
        NOTE: Requires a fortran compiler.
        Parameters:
            met:    String with the type of model to convert to ARL format
                'cmp4': Coamps4 model based on cmp3narl.f90
                            Expects YYYYMMDDHH.tar.gz files in met_dir
            vert_sys:   The vertical coordinate system to use in ARL met files
                        1-pressure sigma
                        2-pressure absolute
                        3-terrain sigma
                        4-hybrid sigma
            met_dir:    Directory with met files to convert
        """
        raise Exception('COAMPS4 preprocessor not included in this version')
        # self.Met = hm.Met(self, met, vert_sys, met_dir)


# ---Functions---
def q_mgr(HyS, HyQ):
    """Manager for the Run Queue that HyRun Processes use.
    Opens the Queue.txt file and grabs and removes the first Run and puts
    in the Run Queue.
    """
    q_file = os.path.join(HyS.var_dir, 'Queue.txt')
    while True:
        run = ''
        lines = open(q_file, 'r').readlines()
        if lines:
            run = lines[0].rstrip()
            open(q_file, 'w').writelines(lines[1:])
        if run:
            HyQ.put(run)
        else:
            break
    return True


# ---Run program main loop--- #FIXME: needs to be moved to __init__ I think
if __name__ == '__main__':
    # If HYSPLITm.py is run normally, running the model only is assumed.
    # Plotting and met processing require importing HYPLITm or running
    # in interactive mode.
    multiprocessing.freeze_support()
    # Initialize an instance of HYSPLITm.Manager()
    Manager = Manager()
    if len(sys.argv) > 1:
        # System flags:
        #   -r: run only, don't add runs
        #   -o: override, rerun any existing requested runs
        if '-r' in sys.argv:
            Manager.start()
        else:
            Manager.hysplit()
    else:
        Manager.hysplit()
