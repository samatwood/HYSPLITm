# plot.py

# HYSPLITm - HYSPLIT Manager

# The HYSPLIT model is maintained by the NOAA Air Resources Lab. Users of this
# program should properly credit and reference NOAA ARL.
# For more information please visit:
#   http://www.arl.noaa.gov/HYSPLIT_info.php
#   http://www.arl.noaa.gov/disclaimer.php

# See the README for more information.

import os
import datetime as dt
import numpy as np
import pylab as pl
import HYSPLITm.model as model
import matplotlib.dates as mdates

# Basemap not always installed
try:
    from mpl_toolkits.basemap import Basemap

    bmf = True
except ImportError:
    bmf = False


class Plot(model.HyBase):
    """Creates various plots of HYSPLIT output.
    # TODO: lots of updates needed in these methods.
    """

    def __init__(self, parent):
        self._parent = parent
        # Setup directories
        self.h_dir = self._parent._h_dir
        self.var_dir = self._parent._var_dir
        self.work_dir = self._parent._work_dir
        # Print warning if Basemap wasn't found
        if not bmf:
            print 'WARNING: Basemap not found in mpl_toolkits. ' + \
                  'Plotting maps will not work.'

    def map_setup(self, map_extents, figsize=(14, 8), r=None,
                  cax=None, zorder=0, **kwargs):
        """Sets up a BaseMap figure to plot trajectories onto.
        Returns the figure or axes object, and map object as tuple, (fig, m)
        Can setup a map an existing figure with cax.
        Parameters:
            map_extents:    Map corners as [ll_lat,ll_lon, ur_lat,ur_lon]
                            Lat(-90:90), Lon(-180:180)
            r:              An example run if no extents given
            figsize:        Figure size
            cax:            The axes to draw the basemap instance on
            zorder:         Will draw map constituents at this zorder
        """
        if map_extents is None and r is None:
            raise Exception('Either map_extents or example run in r is needed')
        if map_extents is None:
            map_extents = self.auto_extents(r, **kwargs)
        if cax is None:
            fig = pl.figure(figsize=figsize)
        else:
            pl.sca(cax)
        m = Basemap(projection='cyl', resolution='l',
                    llcrnrlat=map_extents[0], urcrnrlat=map_extents[2],
                    llcrnrlon=map_extents[1], urcrnrlon=map_extents[3])

        m.drawcoastlines(linewidth=0.5, zorder=zorder)
        m.drawcountries(linewidth=0.5, zorder=zorder)
        m.drawstates(linewidth=0.5, zorder=zorder)
        if cax is None:
            return fig, m
        return cax, m

    def auto_extents(self, r, lat_space=15, lon_space=25, **kwargs):
        """Returns map extents around the receptor for a run name, r.
        """
        var = self.get_var(r)
        lat = var['lat']
        lon = var['lon']
        bottom = lat - lat_space
        top = lat + lat_space
        left = lon - lon_space
        right = lon + lon_space
        if bottom < -90.:
            bottom = -90.
        if top > 90.:
            top = 90.
        if left < -180.:
            left += 360.
        if right > 180:
            right -= 360.
        return [bottom, left, top, right]

    def receptor_grid(self, x, y, grid_deg, xlim=None, ylim=None, recep_center=True):
        """Returns a grid around a receptor.
        x and y are the receptor location in degrees with spacing grid_deg
        and optional max and min limits (each as a list of [min,max]).
        If recep_centr is True, will have the receptor at the center of a grid
        box. If False, the receptor will be at the intesection of 4 grid boxes.
        The grid box edges are returned as a tuple (xedges, yedges).
        """
        if xlim is None:
            xlim = [-180., 180.]
        if ylim is None:
            ylim = [-90., 90.]
        if recep_center:
            nxp = np.ceil((xlim[1] - x + grid_deg / 2.) / grid_deg)
            nxn = np.ceil((x - xlim[0] - grid_deg / 2.) / grid_deg)
            nyp = np.ceil((ylim[1] - y + grid_deg / 2.) / grid_deg)
            nyn = np.ceil((y - ylim[0] - grid_deg / 2.) / grid_deg)

            nx = nxp + nxn + 1
            ny = nyp + nyn + 1

            xe = np.linspace(x - grid_deg / 2. - (grid_deg * nxn),
                             x - grid_deg / 2. + (grid_deg * nxp), nx)
            ye = np.linspace(y - grid_deg / 2. - (grid_deg * nyn),
                             y - grid_deg / 2. + (grid_deg * nyp), ny)
        else:
            nxp = np.ceil((xlim[1] - x) / grid_deg)
            nxn = np.ceil((x - xlim[0]) / grid_deg)
            nyp = np.ceil((ylim[1] - y) / grid_deg)
            nyn = np.ceil((y - ylim[0]) / grid_deg)

            nx = nxp + nxn + 1
            ny = nyp + nyn + 1

            xe = np.linspace(x - (grid_deg * nxn), x + (grid_deg * nxp), nx)
            ye = np.linspace(y - (grid_deg * nyn), y + (grid_deg * nyp), ny)

        if xe[-1] > 180.:
            xe = xe[:-1]
            nx -= 1
        if xe[0] < -180.:
            xe = xe[1:]
            nx -= 1

        if ye[-1] > 180.:
            ye = ye[:-1]
            ny -= 1
        if ye[0] < -180.:
            ye = ye[1:]
            ny -= 1

        return xe, ye

    def traj_plot(self, m, runs,
                  color_type='date',
                  cmap='jet', cb=True, day_int=1,
                  zorder=10,
                  tsmin=None, tsmax=None,
                  scaleminheight=0., scalemaxheight=5000.,
                  scalemaxage=None,
                  daily_ticks=False, tick_m='x', tick_ms=30,
                  hours=None,
                  run_missing=False,
                  fig_type='png', **kwargs):
        """#TODO:
        if cb, will plot colorbar. if -1, will return the scalar mappable for
        creating a colorbar
        zorder default is 10. will draw trajectories at high zorder by default
        """
        # Number of runs to plot
        n_traj = len(runs)
        # Check if each run exists and run if not
        if run_missing:
            q_file = open(os.path.join(self.var_dir, 'Queue.txt'), 'a')
        for name in runs:
            if self.check_ex(name) is False:
                raise Exception
                if run_missing:
                    q_file.write(name + '\n')
        if run_missing:
            if len(open(q_file, 'r').readlines()) > 0:
                os.system('python %s -r' % os.path.join(self.h_dir, 'HYSPLITm.py'))
            q_file.close()

        # Get needed plot parameters from arguments or from run names
        first_run_vars = self.get_var(runs[0])
        cmap = pl.get_cmap(cmap)
        if hours is None:
            hours = np.abs(first_run_vars['dur'])
        if daily_ticks:
            d_ticks = (np.array(range(np.abs(hours) / 24)) + 1) * 24
        else:
            d_ticks = []

        # Get data for all runs
        data, var = self.get_hy_runs(runs, self.HyS.out_dir, metvars=[])
        # trajectory timestamp for date runs
        if color_type == 'date':  # Need to add 2000 to the year
            ts = [dt.datetime(int(data[i, 0, 0]) + 2000, int(data[i, 0, 1]),
                              int(data[i, 0, 2]), int(data[i, 0, 3]),
                              int(data[i, 0, 4])) for i in range(n_traj)]
            ts_ = mdates.date2num(ts)
            if tsmin is None:
                tsmin = min(ts_)
            if tsmax is None:
                tsmax = max(ts_)
            if isinstance(tsmin, dt.datetime):
                tsmin = mdates.date2num(tsmin)
            if isinstance(tsmax, dt.datetime):
                tsmax = mdates.date2num(tsmax)
        # Plot runs
        #   data indicies:
        #       7: lat    8: lon    9: height   6: age
        for i in range(n_traj):
            # plot trajectory, colored by color_type
            if color_type == 'num':
                c_ = cmap(i / (n_traj - 1.))
                self.tplot = m.plot(data[i, :(hours + 1), 8], data[i, :(hours + 1), 7], color=c_, zorder=zorder,
                                    **kwargs)
                for hr in d_ticks:
                    m.scatter(data[i, hr, 8], data[i, hr, 7], tick_ms, marker=tick_m, color=c_, zorder=zorder)
            if color_type == 'date':
                c_ = cmap((ts_[i] - tsmin) / (tsmax - tsmin))
                self.tplot = m.plot(data[i, :(hours + 1), 8], data[i, :(hours + 1), 7], color=c_, zorder=zorder,
                                    **kwargs)
                for hr in d_ticks:
                    m.scatter(data[i, hr, 8], data[i, hr, 7], tick_ms, marker=tick_m, color=c_, zorder=zorder)
            if color_type == 'alt':
                for j in range(hours):
                    c_ = cmap((data[i, j, 9] - scaleminheight) / (scalemaxheight - scaleminheight))
                    self.tplot = m.plot([data[i, j, 8], data[i, j + 1, 8]], [data[i, j, 7], data[i, j + 1, 7]],
                                        color=c_, zorder=zorder, **kwargs)
                    if j in d_ticks:
                        m.scatter(data[i, j, 8], data[i, j, 7], tick_ms, marker=tick_m, color=c_, zorder=zorder)
            if color_type == 'age':
                for j in range(hours):
                    self.tplot = m.plot([data[i, j, 8], data[i, j + 1, 8]], [data[i, j, 7], data[i, j + 1, 7]],
                                        color=cmap(data[i, j, 6] / scalemaxage), zorder=zorder, **kwargs)
                    if j in d_ticks:
                        m.scatter(data[i, j, 8], data[i, j, 7], tick_ms, marker=tick_m, color=c_, zorder=zorder)

        # Plot Color bar #TODO: finish other types
        if cb:
            if color_type == 'date':
                sm = pl.cm.ScalarMappable(norm=pl.Normalize(vmin=tsmin, vmax=tsmax), cmap=cmap)
                sm._A = []
                if cb == -1:
                    return sm
                cb = pl.colorbar(sm, ticks=mdates.DayLocator(interval=day_int), format=mdates.DateFormatter('%d %b'))
                cb.set_label('Date', size=14)
                # pl.title('HYSPLIT %i hr forward trajectories - %i m %s'%(hours,run_height,alt_type_dict[alt_type]))
            if color_type == 'alt':
                sm = pl.cm.ScalarMappable(norm=pl.Normalize(vmin=scaleminheight, vmax=scalemaxheight), cmap=cmap)
                sm._A = []
                if cb == -1:
                    return sm
                cb = pl.colorbar(sm)
                cb.set_label('Height (m)', size=14)
                #pl.title('HYSPLIT %i hr %s %strajectories - %i m %s'%\
                #         (hours,ver_dict[ver],dur_type,run_height,alt_type_dict[alt_type]))

    def rta_plot(self, m, runs,
                 max_endpt=None,
                 grid_deg=0.5,
                 levs=None, clevs=None,
                 nlevs=20, norm_max=True,
                 cmap='jet', cb=True,
                 zorder=-1, ms=40,
                 **kwargs):
        """Contoured Residence Time Analysis plot.
        #TODO: fill in all this stuff fully.
        Parameters:
            m:          Basemap instance
            runs:       A runset to be plotted on m
            max_endpt:  The maximum index to include for each run
                        This is effectively the maximum trajectory age to use
                        assuming endpoints are separated by one hour.
            grid_deg:   The size of grid box edges to create
            levs:       Fixed levels to use.
                        If not None, will use these and override nlevs.
            norm_max:   Normalized grid box frequencies by fraction of maximum
                        residence time if True. Otherwise, will normalize as
                        fraction of total residence time (i.e. num endpoints).

        """
        # Check if each run exists
        for name in runs:
            if self.check_ex(name) is False:
                raise Exception('Run not found: %s' % name)

        # Get data for all runs
        data, var = self.get_hy_runs(runs, self.HyS.out_dir, metvars=[])
        ntraj, nendpt = data.shape[0:2]
        if max_endpt is None:
            max_endpt = nendpt
        else:
            if max_endpt > nendpt - 1:
                raise Exception('Not enough enpoints for this max_endpt')
        # Get lat and lon for each endpoint
        latloc = var.index('lat')
        lonloc = var.index('lon')
        e_lat = []
        e_lon = []
        for nt in range(ntraj):
            for ne in range(1, max_endpt):
                e_lat.append(data[nt, ne, latloc])
                e_lon.append(data[nt, ne, lonloc])
        e_lat = np.array(e_lat)
        e_lon = np.array(e_lon)

        # Create grid around receptor
        r_x = data[0, 0, lonloc]
        r_y = data[0, 0, latloc]
        xedges, yedges = self.receptor_grid(data[0, 0, lonloc], data[0, 0, latloc],
                                            grid_deg, recep_center=True,
                                            xlim=[e_lon.min(), e_lon.max()],
                                            ylim=[e_lat.min(), e_lat.max()])
        # Create histogram of lat and lon endpoints
        h, xe_, ye_ = np.histogram2d(e_lon, e_lat, bins=[xedges, yedges])
        # Grid centers
        x = (xedges[1:] - xedges[:-1]) / 2. + xedges[:-1]
        y = (yedges[1:] - yedges[:-1]) / 2. + yedges[:-1]
        # Normalize grid box frequencies by:
        if norm_max:
            # fraction of maximum residence time
            rta = h.transpose() / h.max()
            extend = 'neither'
        else:
            # Normalize the frequencies
            rta = h.transpose() / float(e_lat.size)
            extend = 'max'
        if levs is None:
            levs = np.linspace(rta.min(), rta.max(), nlevs + 1)[1:]
        grid = pl.meshgrid(x, y)
        m.scatter([r_x], [r_y], s=ms, c='k', marker='x', linewidth=2)
        m.cplot = m.contourf(grid[0], grid[1], rta, levels=levs, colors=clevs,
                             extend=extend, cmap=cmap, zorder=zorder)
        if cb:
            m.colorbar(m.cplot, size='3%')
