#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

import sys
sys.path.append('.')

import numpy as np
import logging
import cPickle
import matplotlib as mpl
from scipy.io import readsav

from load_frames import load_mdsframes
from misc.phantom_helper import make_rz_array, find_sol_pixels
from blob_tracking import blob_tracking as blob_tracking_fun
from plotting import plot_trail_geom

from analysis.velocity_analysis import velocity_analysis


frame0 = 0
nframes = 10000
save = False

# Directory, where frames from mds are stored in
#datadir = '/Volumes/Backup/cmod_data/phantom'
datadir = '/Users/ralph/uni/cmod/tmp_data/'
# Directory where log files and plots are written to
outdir = '/Users/ralph/uni/cmod'
logfile = '%s/blob_tracking/blobtracking.log' % (outdir)

# Store information about density scan with loggin facility
logger = 0.
logger = logging.getLogger('mylogger')
logger.setLevel(logging.INFO)
file_handle = logging.FileHandler(logfile)
# file_handle.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
# file_handle.setFormatter(formatter)
logger.addHandler(file_handle)

logger.info('Start logging!!!')

shotlist = [1111208020]
thresh = 2.5

# Identify peaks in this amplitude range
minmax = [thresh, 10.0]
# Identify peaks in this trigger box
trigger = np.array([40, 50, 16, 48])

for shotnr in shotlist:
    log_str = '*********************** Shot #%10d ' % shotnr
    try:
        logger.info(log_str)
    except:
        print log_str
    # Load frames
    frames, fi = load_mdsframes(shotnr, path=datadir, fname='%10d_testframes2.npz' % shotnr, varname='frames_normalized')
    # Load separatrix data
    fname_sep = '%s/%10d_separatrix.sav' % (datadir, shotnr)
    sep = readsav(fname_sep)

    sol_px = find_sol_pixels(sep)

    print type(sol_px)
    print sol_px.shape
    print sol_px
    
    raise ValueError
    # #########################################################################
    # Do blob tracking
    #
    trails = blob_tracking_fun(shotnr, frames, trigger, minmax)

    print 'Found %d blob trails' % len(trails)

    # #########################################################################
    # Compute geometry data from GPI and plot some blob trails
    #
    rz_array, transform_data = make_rz_array(fi)
    xxi, yyi = np.meshgrid(np.linspace(rz_array[:, :, 0].min(), rz_array[:, :, 0].max(), 64),
                           np.linspace(rz_array[:, :, 1].min(), rz_array[:, :, 1].max(), 64))
    xyi = np.concatenate((xxi[:, :, np.newaxis], yyi[:, :, np.newaxis]), axis=2)
    # Compute the indices between LCFS and limiter shadow

    #domain = np.array(good_idx)[:, :].tolist()

    #for trail in trails:
    #    # Plot a simple trail
    #    #plot_trail_simple(trail, frames, plot_shape=True)
    #    # Plot trail with GPI geometry overlay
    #    plot_trail_geom(trail, frames, rz_array=rz_array, xyi=xyi, 
    #            trigger_box=trigger, sep_data=sep, 
    #            plot_com=True, plot_max=False, plot_geom=True) 

    # #########################################################################
    # Save trails for later use
    #
    
    #picklefile = open('%d/%d_trails_thresh%2d.pkl' % (shotnr, shotnr,
    #                                                  int(10 * thresh)),
    #                  'wb')
    #cPickle.dump(trails, picklefile, -1)
    #picklefile.close()

    # #########################################################################
    # Compute blobtrail statistics
    #
    velocity_analysis(trails, frames, sol_px, rz_array, xyi)

    #print 'Computing trail statistics'
    #blob_amps, blob_ell, blob_vel, blob_shape, blobs_used_good_domain, fail_list = blob_statistics.statistics_blob_sol(shotnr, trails, fi, frames=frames, good_domain=domain, logger=logger)

    #blob_statistics.plot_blob_stat2(blob_amps, blob_ell, blob_vel, blob_shape,
    #                                fail_list, shotnr, len(trails),
    #                                len(fail_list), save=False, logger=logger)

    #blob_statistics.plot_blob_scatter1(blob_amps, blob_ell, blob_vel,
    #                                   blob_shape, shotnr, fail_list,
    #                                   save=False, logger=logger)

    # ########################################################################
    # 
    #

#    np.savez('shot_100819/ds_thresh%2d_results_20120709.npz' % (int(10*thresh)), blob_amps = blob_amps, blob_ell = blob_ell, blob_vel = blob_vel, blob_shape = blob_shape, fail_list = fail_list )
