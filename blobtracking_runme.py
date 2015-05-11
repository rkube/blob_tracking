#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

import numpy as np
import logging
import cPickle
import matplotlib as mpl
#import itertools
#mpl.use('AGG')
import matplotlib.pyplot as plt
from misc.load_mdsframes import load_mdsframes
from misc import helper_functions
from blob_tracking.blob_tracking import blob_tracking as blob_tracking_fun
from blob_tracking import blob_statistics

frame0 = 0
save = False

# Store information about density scan with loggin facility
logger = 0.
logger = logging.getLogger('density_scan_120711')
logger.setLevel(logging.INFO)
file_handle = logging.FileHandler('logs/test666.log')
# file_handle.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
# file_handle.setFormatter(formatter)
logger.addHandler(file_handle)

logger.info('Start logging!!!')

shotlist = [1100803020]
thresh = 2.5

for shotnr in shotlist:
    print 'Processing shot# %d' % shotnr
    logger.info('Processing shot %d' % shotnr)

    frames, fi = load_mdsframes(shotnr, test=False)
    frames = frames[frame0:, :, :]
    frames[:, 0, 0] = frames.min()
    frames[:, 1, 0] = frames.max()

    # Compute the indices between LCFS and limiter shadow
    good_idx = helper_functions.find_sol_pixels(shotnr, fi)
    domain = np.array(good_idx)[:, :].tolist()

    try:
        print 'Loading detected trails from file'
        picklefile = open('%d/%d_trails_FU_thresh%2d.pkl' %
                          (shotnr, shotnr, int(10 * thresh)), 'rb')
        trails = cPickle.load(picklefile)
        picklefile.close()

        print 'Unpickled blob trails'
        logger.info('Loaded %d blob trails from file' % len(trails))

    except:
        print 'Running blob detection.'
        trails = blob_tracking_fun(shotnr, frames, fi,
                                   minmax=[thresh, 10.0], logger=logger)
        logger.info('Blob detection found %d blobs' % len(trails))
        print 'Found %d blobs' % len(trails)
        picklefile = open('%d/%d_trails_thresh%2d.pkl' % (shotnr, shotnr,
                                                          int(10 * thresh)),
                          'wb')
        cPickle.dump(trails, picklefile, -1)
        picklefile.close()

    print 'Computing trail statistics'
    blob_amps, blob_ell, blob_vel, blob_shape, blobs_used_good_domain, fail_list = blob_statistics.statistics_blob_sol(shotnr, trails, fi, frames=frames, good_domain=domain, logger=logger)

    blob_statistics.plot_blob_stat2(blob_amps, blob_ell, blob_vel, blob_shape,
                                    fail_list, shotnr, len(trails),
                                    len(fail_list), save=False, logger=logger)

    blob_statistics.plot_blob_scatter1(blob_amps, blob_ell, blob_vel,
                                       blob_shape, shotnr, fail_list,
                                       save=False, logger=logger)

#    np.savez('shot_100819/ds_thresh%2d_results_20120709.npz' % (int(10*thresh)), blob_amps = blob_amps, blob_ell = blob_ell, blob_vel = blob_vel, blob_shape = blob_shape, fail_list = fail_list )
