#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

"""
============
tracker
============


Implements blob tracking.
The idea is to identify one coheren large amplitude structure in an image.
Using a flooding and labelin algorithm from pymorph, the coherent region
is identified. Identify coherent regions of similar amplitude in successive frames
and track the motion of this region over a maximum number of frames, or
until no sufficiently large region can be identified


"""

import numpy as np
from geometry import com
from helper_functions import find_closest_region, width_gaussian, TrackingError

def tracker_geom(frames, x0, event, thresh_amp, thresh_dist, blob_ext,
                 plots=False, verbose=False):
    """
    Track the blob in a dynamic intervall forward or backward in time, as long
    as its amplitude is over a given threshold and the peak has detected less
    than a given threshold over consecutive frames.
    The maximum number of frames the blob is tracked for, is given by dim0 of
    frames

    Input:
        frames:     The frames on which we track the centroid motion
        x0:         Original position of the centroid
        event:      ndarray, [I0, t0, R0, z0] Index of original feature to
                    track
        thresh_amp: Threshold for amplitude decay relative to frame0
        thresh_dist Threshold for blob movement relative to previous frame
        blob_ext:   Extend of the blob used for determining its average shape

    Returns:
        numframes:      Number of frames the blob was tracked
        xycom:          COM position of the blob in each frame
        amp:            Amplitude at COM in each frame
        fwhm_rad_idx:   Indices that mark left and right FWHM of the blob
        fwhm_pol_idx:   Indices that mark the lower and upper FWHM of the blob
        blob:           Array that stores the blob extend
    """

    assert (blob_ext % 2 == 0)

    if (verbose is True):
        print 'Called tracker with '
        print '\tevent = ', event
        print '\tthresh_amp = ', thresh_amp
        print '\tthresh_dist = ', thresh_dist
        print '\tblob_ext = ', blob_ext
        print '\tplots = ', plots

    # Maximum number of frames the blob is tracked for is given by
    # dimension 0 of frames
    # tau_max = np.shape(frames)[0]
    tau_max = frames.shape[0]
    I0_last, z0_last, R0_last = event[0], x0[0], x0[1]

    print 'tau_max = %d' % (tau_max)

    # I0 is the threshold amplitude we use for detecting blobs
    # i.e. for blob tracking, we identify later all connected regions that are larger than 
    # I0 * thresh_amp

    if (verbose):
        verb_msg = 'Tracking blob, x = %d, y = %d, I0 = %f' % (R0_last, z0_last, I0_last)
        verb_msg += ' thresh_amp = %f, thresh_dist = %f' % (thresh_amp * I0_last, thresh_dist)
        print verb_msg

    xycom = np.zeros([tau_max, 2])                      # Return values: COM position of peak
    xymax = np.zeros([tau_max, 2])                      # Position of the blob peak
    width_pol = np.zeros(tau_max, dtype='float64')      # Width of Gaussian fit, poloidal direction  
    width_rad = np.zeros(tau_max, dtype='float64')      # Width of Gaussian fit, radial direction
    fwhm_pol = np.zeros([tau_max, 2], dtype='int')      # Poloidal FWHM
    fwhm_rad = np.zeros([tau_max, 2], dtype='int')      # Radial FWHM
    amp = np.zeros([tau_max])                           # Amplitude at COM position

    good_blob = True
    tau = 0     # The offset, where we try to find a good blob
    while (good_blob and tau < tau_max - 1):
        tau += 1                              # Advance indices
        # Load next frame
        event_frame = frames[tau, :, :]
        
        if (verbose):
            print ''
            print 'tau: %d, Last blob position: x=%d, y=%d, I0=%f' %\
                (tau, R0_last, z0_last, I0_last)
    
        # Find the closest centroid region
        try:
            x0, event_masked = find_closest_region(event_frame, thresh_amp * I0_last, 
                    [z0_last, R0_last], max_dist=thresh_dist, verbose=False)
        except TrackingError:
            print 'tau = %d: Lost blob track' % tau
            tau -= 1
            break

        if (verbose):
            print '        New blob position: x=%d, y=%d' % (x0[0], x0[1])

        # Compute maximum and COM position on the new centroid
        xymax[tau, :] = np.unravel_index(event_masked.argmax(), event_frame.shape)
                                         
        # When used to index frames[:,:,:]:
        #     xycom[tau,:] = [index for axis 1, index for axis 2]
        # To be consistent with indexing from xymax, flip this array
        # COM returns com along second dimension at index 0
        xycom[tau, ::-1] = com(event_masked)
        ycom_off, xcom_off = xycom[tau, :].round().astype('int')

        # Follow the peak
        z0_last, R0_last = xymax[tau, :].astype('int')
        amp[tau] = event_frame[xymax[tau, 0], xymax[tau, 1]]

        if verbose:
            vrb_msg = '          New blob position:'
            vrb_msg += '          x_max = (%d,%d)' %(xymax[tau, 0], xymax[tau, 1])
            vrb_msg += '          x_com = (%d,%d)' %(xycom[tau, 0], xycom[tau, 1])
            print vrb_msg

        # Fit a Gaussian on the cross section around the maximum
        sigma, dummy1, dummy2 = width_gaussian(event_frame[xymax[tau, 0], :], xymax[tau, 1])
        width_pol[tau] = sigma
        sigma, dummy1, dummy2 = width_gaussian(event_frame[:, xymax[tau, 1]], xymax[tau, 0])
        width_rad[tau] = sigma

        if (plots):
            plt.figure(figsize=(18, 6))
            plt.subplot(131)
            plt.title('frame %d' % (tau))
            plt.contourf(event_frame, 64, cmap=plt.cm.hot)
            plt.plot(xycom[tau, 1], xycom[tau, 0], 'ko')
            plt.plot(xymax[tau, 1], xymax[tau, 0], 'k^')
            plt.colorbar()
            plt.xlabel('x / px')
            plt.ylabel('y / px')

            plt.subplot(132)
            plt.plot(event_frame[xycom[tau, 0], :], 'b', label='ax0')
            plt.plot(xycom[tau, 0], event_frame[xycom[tau, 0], xycom[tau, 1]], 'ko')
            plt.plot(event_frame[:, xycom[tau, 1]], 'g', label='ax1')
            plt.plot(xycom[tau, 1], event_frame[xycom[tau, 0], xycom[tau, 1]], 'ko')
            plt.title('COM crosssections')

            plt.subplot(133)
            plt.plot(event_frame[xymax[tau, 0], :], 'b', label='ax0')
            plt.plot(xymax[tau, 0], event_frame[xymax[tau, 0], xymax[tau, 1]], 'k^')
            plt.plot(event_frame[:, xymax[tau, 1]], 'g', label='ax1')
            plt.plot(xymax[tau, 1], event_frame[xymax[tau, 0], xymax[tau, 1]], 'k^')
            plt.title('MAX crosssections')

    # end while (good_blob and tau < tau_max):

    # Crop tracking results to actual length
    amp = amp[1:tau + 1]
    xycom = xycom[1:tau + 1, :]
    xymax = xymax[1:tau + 1, :]
    width_rad = width_rad[1:tau + 1]
    width_pol = width_pol[1:tau + 1]

    if (plots):
        plt.show()

    return tau, amp, xycom, xymax, width_rad, width_pol



# End of file tracker.py
