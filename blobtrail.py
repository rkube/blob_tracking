#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

"""
=========
blobtrail
=========

.. codeauthor :: Ralph Kube <ralphkube@gmail.com>

A class that defines a blob event in a sequence of frames 
from 2d turbulence imaging
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import griddata
from scipy.optimize import leastsq
from tracker import tracker_geom
from geometry import com 
from helper_functions import find_closest_region, width_gaussian


class blobtrail:
    """
    A realization of a blob event and its associated trail
    input:
    =====

    * frames:      imaging frames. axis0: time
                                   axis1: poloidal direction
                                   axis2: radial direction
    * event:       structured array (intensity, t0, z0, R0)
    * shotnr:      shotnr
    * tau_max:     Maximal number of frames for fwd/bwd tracking
    * thresh_dist: Max. distance a peak is allowd to travel per frame
    * fwhm_max_idx:Maximum width allowed, in pixel
    * blob_ext:    Size of the blob
    * thresh_amp:  Threshold in percentage of original amplitude before a
                   blob is considered to be lost by the tracker
    * frame0:      offset frame for frames array
    * doplots:     Plot output
    """

    def __init__(self, frames, event, shotnr, tau_max=7,
                 thresh_dist=8., fwhm_max_idx=10, blob_ext=8, thresh_amp=0.6,
                 frame0=0, doplots=False):

        # Store all the parameters this blob has been tracked with.
        self.event = event                  #
        self.tau_max = tau_max              #
        self.thresh_dist = thresh_dist      #
        self.fwhm_max_idx = fwhm_max_idx    #
        self.blob_ext = blob_ext            #
        self.thresh_amp = thresh_amp        #
        self.frame0 = frame0                #
        self.dt = 2.5e-6                    #
        self.shotnr = shotnr                #

        # Error flags that signal something went wrong with blob tracking
        self.invalid_fw_tracking = False
        self.invalid_bw_tracking = False

        # Number of frames the blob is tracked backwards. Updated in track_backward
        self.tau_b = 0
        # Number of frames the blob is tracked forwards. Updated in track_forward
        self.tau_f = 0
        # Find the first centroid position
        frame0 = frames[tau_max, :, :]
        x0, dummy = find_closest_region(frame0, thresh_amp * event[0], 
                                        [event[2], event[3]], 
                                        max_dist = 5.0, verbose=False)
        # Compute maximum and COM coordinates of the blob within in the first centroid
        xycom0 = com(dummy)
        xycom0 = np.array(xycom0)
        #print 'before: ', xycom0
        # See com
        #xycom0 = xycom0[::-1]
        #print 'after', xycom0
        
        xymax0 = np.unravel_index(dummy.argmax(), frames[0, :, :].shape)
        xymax0 = np.array(xymax0)

        sigma_rad0, dummy1, dummy2 = width_gaussian(frame0[xymax0[0], :], xymax0[1])
        sigma_pol0, dummy1, dummy2 = width_gaussian(frame0[:, xymax0[1]], xymax0[0])

        # Track blob forwards and backwards, combine results
        self.track_forward(frames, x0, doplots=doplots)
        self.track_backward(frames, x0, doplots=doplots)

        # Combine results from forward and backward tracking
        # Values about the blob path [backward, 0, forward]
        self.tau = np.arange(-self.tau_b, self.tau_f + 1)
        self.amp = np.concatenate([self.amp_b[::-1], [event[0]], self.amp_f])
        self.xycom = np.concatenate([self.xycom_b[::-1], xycom0[np.newaxis, :], self.xycom_f])
        self.xymax = np.concatenate([self.xymax_b[::-1], xymax0[np.newaxis, :], self.xymax_f])
        self.ell_rad = np.concatenate([self.width_rad_b[::1], np.array([sigma_rad0]), self.width_rad_f])
        self.ell_pol = np.concatenate([self.width_pol_b[::-1], np.array([sigma_pol0]), self.width_pol_f])

    def track_backward(self, frames, x0, doplots=False):
        """
        Track blob backward frame0 to beginning of frames
        """
        # Reverse the frame order, so that forward tracking is tracking backwards
        # i.e.
        # frames[-tau_max ... tau_max, :, :]
        # pass the time index to tracker, so that it has
        # frames[[0, -1, -2, .., -tau_max], :, :]
        idx_frames = self.tau_max + np.arange(0, -self.tau_max - 1, -1, dtype='int')

        res = tracker_geom(frames[idx_frames, :, :], x0, self.event,
                           self.thresh_amp, self.thresh_dist, self.blob_ext,
                           plots=doplots, verbose=True)
        self.tau_b, self.amp_b, self.xycom_b, self.xymax_b, self.width_rad_b, self.width_pol_b = res


    def track_forward(self, frames, x0, doplots=False):
        """
        Track blob forward from frame0
        """
        #print 'Tracking forward...'
        res = tracker_geom(frames[self.tau_max:, :, :], x0, self.event,
                           self.thresh_amp, self.thresh_dist, self.blob_ext,
                           plots=doplots, verbose=False)
        self.tau_f, self.amp_f, self.xycom_f, self.xymax_f, self.width_rad_f, self.width_pol_f = res


    def get_frame0(self):
        """
        The index for the frame where the blob was detected
        """
        return self.frame0

    
    def length(self):
        """
        Length of the trail
        """
        return self.tau.size


    def get_xycom(self):
        """
        Return the COM position of the blob
        """
        return self.xycom

    def get_xymax(self):
        """
        Return the MAX position of the blob
        """
        return self.xymax


    def get_fwhm_pol(self):
        """
        Return the previously computed poloidal width of the blob
        """
        return self.fwhm_ell_pol


    def get_fwhm_rad(self):
        """
        Return the previously computed radial width of the blob
        """
        return self.fwhm_ell_rad


    def get_ell_pol(self):
        """
        Return error from length fitting
        """
        return self.ell_pol

    def get_ell_rad(self):
        """
        Return error from length fitting
        """
        return self.ell_rad


    def get_amp(self):
        """
        Return the amplitude (maximum intensity) of the blob
        """
        return self.amp


    def get_tau(self):
        """
        Return the frames in blob trail relative to the frame number where the blob was detected
        """
        return self.tau

    def get_event(self):
        """
        Return the event where we started tracking from
        """
        return self.event

    def get_event_frames(self):
        """
        Return the frames in which the blob event occurs
        """
        return self.tau + self.event[1]

    def get_blob_shape(self, frames, frameno = None, position = 'COM'):
        """
        Return a the shape of the blob centered around its COM position
        position:   Return blob position at center of mass ('COM') or maximum ('MAX')
        frameno:    Returns the blob shape in the specified range, this range must be within [-tau_b : tau_f]
        """

        assert( position in ['COM', 'MAX'] )

        if (frameno is not None) :
            assert (isinstance(frameno, np.ndarray))
            assert (frameno.max() <= self.tau_f)
            assert (frameno.min() >= -self.tau_b)

            blob_shape = np.zeros([np.size(frameno), 2 * self.blob_ext, 2 * self.blob_ext])
            t_off = frameno

        else:
            blob_shape = np.zeros([np.size(self.tau), 2 * self.blob_ext, 2 * self.blob_ext])
            t_off = np.arange(np.size(self.tau))


        if (position is 'COM'):
            x_off, y_off = self.xycom[:,0].astype('int'), self.xycom[:,1].astype('int')
        elif (position is 'MAX'):
            x_off, y_off = self.xymax[:,0].astype('int'), self.xymax[:,1].astype('int')

        for t_idx, t in enumerate(t_off):
            blob_shape[t_idx, :, :] = frames[t + self.event[1] + self.frame0, y_off[t_idx] - self.blob_ext : y_off[t_idx] + self.blob_ext, x_off[t_idx] - self.blob_ext : x_off[t_idx] + self.blob_ext]

        print 'blob_shape finished'
        return blob_shape

# End of file blobtrail.py
