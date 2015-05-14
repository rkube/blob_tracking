#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

import numpy as np
import logging
import cPickle
import matplotlib as mpl
#from misc.load_mdsframes import load_mdsframes
from load_mdsframes import load_mdsframes
import helper_functions
from blob_tracking import blob_tracking as blob_tracking_fun
#from blob_tracking.blob_tracking import blob_tracking as blob_tracking_fun
#from blob_tracking import blob_statistics

frame0 = 10000
nframes = 10000
save = False

# Directory, where frames from mds are stored in
datadir = '/Volumes/Backup/cmod_data/phantom'
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

shotlist = [1100803020]
thresh = 2.5

minmax = [thresh, 10.0]
trigger = np.array([40, 50, 16, 48])

# Create coordinate mapping for the pixels
# This is dummy code, replace with the true coordinate mapping of the input data
xrg = np.arange(64)
xx, yy = np.meshgrid(xrg, xrg)
rz_array = np.zeros([64, 64, 2])
rz_array[:, :, 0] = xx[:, :]
rz_array[:, :, 0] = yy[:, :]

for shotnr in shotlist:
    print 'Processing shot# %d' % shotnr
    logger.info('Processing shot %d' % shotnr)

    # Directory where we find the frames for the current shot
    ddir_shot = '%s/%10d' % (datadir, shotnr)

    frames, fi = load_mdsframes(shotnr, test=False, path=ddir_shot)
    frames = frames[frame0:frame0 + nframes, :, :]
    print 'loaded frames, min = %f, min = %f' % (frames.min(), frames.max())
    # Load rz_array for frames
    print '********* Using dummy rz_array **********'
    #rz_array, transform_data = make_rz_array(fi)
    #frames[:, 0, 0] = frames.min()
    #frames[:, 1, 0] = frames.max()

    # Compute the indices between LCFS and limiter shadow
    good_idx = helper_functions.find_sol_pixels(shotnr, fi, datadir='./data')
    domain = np.array(good_idx)[:, :].tolist()


    print 'Running blob detection.'
    trails = blob_tracking_fun(shotnr, frames, rz_array,
                               trigger, minmax, logger)

    #trails = blob_tracking_fun(shotnr, frames, fi,
    #                           minmax=[thresh, 10.0], logger=logger)
    #logger.info('Blob detection found %d blobs' % len(trails))
    #print 'Found %d blobs' % len(trails)
    #picklefile = open('%d/%d_trails_thresh%2d.pkl' % (shotnr, shotnr,
    #                                                  int(10 * thresh)),
    #                  'wb')
    #cPickle.dump(trails, picklefile, -1)
    #picklefile.close()

    #print 'Computing trail statistics'
    #blob_amps, blob_ell, blob_vel, blob_shape, blobs_used_good_domain, fail_list = blob_statistics.statistics_blob_sol(shotnr, trails, fi, frames=frames, good_domain=domain, logger=logger)

    #blob_statistics.plot_blob_stat2(blob_amps, blob_ell, blob_vel, blob_shape,
    #                                fail_list, shotnr, len(trails),
    #                                len(fail_list), save=False, logger=logger)

    #blob_statistics.plot_blob_scatter1(blob_amps, blob_ell, blob_vel,
    #                                   blob_shape, shotnr, fail_list,
    #                                   save=False, logger=logger)

#    np.savez('shot_100819/ds_thresh%2d_results_20120709.npz' % (int(10*thresh)), blob_amps = blob_amps, blob_ell = blob_ell, blob_vel = blob_vel, blob_shape = blob_shape, fail_list = fail_list )
