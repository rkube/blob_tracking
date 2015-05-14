#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
#from misc.phantom_helper import make_rz_array
# from blob_tracking.detect_peak import detect_peak_3d
import blobtrail
from detect_peak import detect_peak_3d
# from scipy.interpolate import griddata
# import cPickle as pickle
from scipy.io import readsav


"""
Run blob detection on a set of GPI frames and store information about the
blob in a blobtrail object. Return the list of blobtrail objects
"""


def blob_tracking(shotnr, frames, rz_array, trigger=np.array([40, 50, 16, 48]), 
                  minmax=np.array([2.0, 10.0]),
                  logger=None, 
                  datadir='/Volumes/Backup/cmod_data/phantom/',
                  frame0=0):

    """
    Run blob detection and tracking on the provided frames
    shotnr:   integer, 10 digits. Shot identifier 
    frames:   ndarray, axis0: time, axis1: radial, axis2: poloidal
    rz_array: ndarray, Maps pixel space to coorinate space
                       axis0: pixel coordinate, 
                       axis1: pixel coordinate, 
                       axis2[0]: x -coordinate, axis2[1]: y-coordinate
    trigger:  ndarray, Array in which blobs are to be detected: [r_lo, r_up, z_low, z_up]
    minmax:   ndarray, minimum / maximum threshold for blob detection
    logger:   logger facility to use
    datadir:  Optional, where data is located
    frame0:   Optional, frame offset to store in blobtrail object
    """

    np.set_printoptions(linewidth=999999)
    # Peaks within 2.5 and 10.0 times the rms
    minmax = np.array(minmax)
    # Deadtime after a peak in which no blob is detected
    lag = 20
    # Total frame offset used in this script.
    # toffset = frame0 + lag
    # Maximal frames for which we track a blob
    tau_max = 7
    # 1 frame is 2.5µs
    # dt = 2.5e-6

    log_msg = 'Starting blob tracking for shot %d'
    try:
        logger.info(log_msg)
    except:
        print log_msg

    # Load separatrix data for shot
    #s = readsav('%s/%d/%d_separatrix.sav' % (datadir, shotnr, shotnr), verbose=False)

    # Detect peaks
    # The detect_peak_3d returns time indices of blob events relative for
    # the array passed to it. Remember to add frame0 to idx_event[t0] to
    # translate to the frame indexing used in this script.
    idx_events = detect_peak_3d(frames, trigger, minmax, lag, rel_idx=False)
    print idx_events
    print frames.shape
    num_events = np.shape(idx_events)[0]

    log_msg = '%d blob events detected' % (num_events)
    try:
        logger.info(log_msg)
    except:
        print log_msg

    # Define the events we will analyze
#    event_range = np.arange(num_events)
    event_range = np.arange(46, 47)
#    num_events  = np.size(event_range)

    # Get R,z projection, grid data
    #rz_array, transform_data = make_rz_array(frame_info)
    xxi, yyi = np.meshgrid(np.linspace(rz_array[:, :, 0].min(),
                                       rz_array[:, :, 0].max(), 64),
                           np.linspace(rz_array[:, :, 1].min(),
                                       rz_array[:, :, 1].max(), 64))
    xyi = np.concatenate((xxi[:, :, np.newaxis],
                          yyi[:, :, np.newaxis]), axis=2)
    trails = []
    fail_list = []
    failcount = 0

    for idx, event in enumerate(idx_events[event_range]):
        # I0 = event[0]
        t0 = event[1] 
        # z0 = event[2]
        # R0 = event[3]

        print ''
        print '=============================================================='
        print 'Tracking peak %d / %d, frame %d' % (idx, num_events, t0)
        # try:
        #plt.figure()
        #plt.contourf(frames[t0, :, :], 64)
        #plt.colorbar()
        #plt.show()
        #print 'frames.max  = %f' % frames[t0, :, :].max()
        newtrail = blobtrail.blobtrail(frames[t0 - tau_max:
                                              t0 + tau_max + 1, :, :],
                                       event, shotnr,
                                       thresh_amp=0.7, blob_ext=14,
                                       thresh_dist=8.,
                                       fwhm_max_idx=18,
                                       frame0=frame0,
                                       doplots=False)

        print ''
        print 'blob trail max:', newtrail.get_trail_max()
        print 'blob trail COM:', newtrail.get_trail_com()
        print ''
        print '=========================================================='

        if (np.size(newtrail.get_tau()) < 4):
            fail_list.append(idx)
            failcount += 1

            log_str = 'Peak %d: Trail too short: %d frames' %\
                (idx, newtrail.get_tau.size)

            try:
                logger.info(log_str)
            except:
                print log_str

            continue

        #except ValueError, IndexError:
        #    fail_list.append(idx)
        #    failcount += 1
        #    log_str = 'Failed to track blob %d / %d' % (idx, num_events)
        #    try:
        #        logger.info(log_str)
        #    except:
        #        print log_str
        #    continue

        print 'Computing blob width'
        try:
            newtrail.compute_width_gaussian(frames, rz_array, position='MAX',
                                            i_size_max=10, plots=False)

        except:
            fail_list.append(idx)
            failcount += 1
            log_str = 'Peak %d: Unable to compute FWHM' % (idx)
            try:
                logger.info(log_str)
            except:
                print log_str

        newtrail.plot_trail(frames, rz_array=rz_array, xyi=xyi,
                            trigger_box=trigger, sep_data=None,
                            plot_com=True, plot_shape=True, plot_geom=False,
                            save_frames=False)
        trails.append(newtrail)

    return trails

# End of file blob_tracking.py
