#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

"""
================
Helper functions
================

.. codeauthor :: Ralph Kube <ralphkube@gmail.com>
Helper functions to make blob analysis and tracking easier

* com......Computes the center of mass position for a 2d array
* com_Rz...Computes the center of mass position for a 2d array with R and z
           coordinates
* fwhm.....Compute the full width half maximum of a 1d array

"""

import numpy as np
import pymorph as pm
import matplotlib.pyplot as plt
# from phantom_helper import make_rz_array
# from scipy.interpolate import griddata
from scipy.optimize import leastsq
from scipy.io import readsav



def frac_to_idx(frac, f_min, f_max, nbins):
    """
    Given a float within in the interval [f_min,f_max], return the index in
    which of nbins equidistant bins it belongs
    """
    return np.floor((frac - f_min) * float(nbins) /
                    (f_max - f_min)).astype('int')


def com(array, xx=None, yy=None):
    """
    Return the center of mass position of the array x_com and y_com
    x_com = int(x * array) / int(array),
    y_com = int(y * array) / int(array)
    If xx and yy are not specified, use x = 0,1,...,np.shape(array)[0],
    and y = 0, 1, ..., np.shape(array)[1]

    Returns
        x_com, y_com
    """
    array = array.astype('float')
    if (xx is None and yy is None):
        xx, yy = np.meshgrid(np.arange(0, array.shape[0], 1.0),
                             np.arange(0, array.shape[1], 1.0))
        # If dx and dy are not specified, assume a regularly spaced grid
        return ((xx * array).sum() / array.sum(),
                (yy * array).sum() / array.sum())
    else:
        # Compute the increments in x and y
        dx = np.zeros_like(xx)
        dy = np.zeros_like(yy)
        dx[:, :-1] = xx[:, 1:] - xx[:, :-1]
        dx[:, -1] = dx[:, -2]

        dy[:-1, :] = yy[1:, :] - yy[:-1, :]
        dy[-2, :] = dy[-1, :]
        # Surface element
        dA = np.abs(dx) * np.abs(dy)
        return ((xx * array * dA).sum() / (array * dA).sum(),
                (yy * array * dA).sum() / (array * dA))

#    try:
#        dx = (xx[:,1:] - xx[:,:-1]).max()
#        dy = (yy[1:,:] - yy[:-1,:]).max()
#        return ((xx * array).sum() / array.sum(),
#                (yy * array).sum() / array.sum()
#    except TypeError:
#        xx, yy = np.meshgrid(np.arange(0, np.shape(array)[0]),
#                             np.arange(0, np.shape(array)[1]))
#        dx, dy = 1
#        #dx = np.max( xx[:,1:] - xx[:,:-1] )
#        #dy = np.max( yy[1:,:] - yy[:-1,:] )
#        #return (xx * array * dx * dy).sum() / (array * dx * dy).sum(), (yy *
#        array * dx * dy).sum() / (array * dx * dy).sum()


def com_rz(array, RR, zz):
    """
    Return the center of mass position on the irregulary spaced RR, zz array:
    R_com = int ( R * n * dA ) / int ( n * dA ), along second dimension
    z_com = int ( z * n * dA ) / int ( n * dA ), along first dimension
    """
    array = array.astype("float")
    dR, dz = np.zeros_like(RR), np.zeros_like(zz)
    dR[:, :-1] = RR[:, 1:] - RR[:, :-1]
    dR[:, -1] = dR[:, -2]
    dz[:-1, :] = zz[1:, :] - zz[:-1, :]
    dz[-1, :] = dz[:, -2]

    dA = np.abs(dR) * np.abs(dz)

    # COM along second dimension, COM along first dimension
    return((RR * array * dA).sum() / (array * dA).sum(),
           (zz * array * dA).sum() / (array * dz).sum())
    # return np.sum( RR * array * dR * dz ) / np.sum (array * dR * dz ) ,\
    #     np.sum( zz * array * dR * dz ) / np.sum (array * dR * dz )


def fwhm(array):
    """
    Computes the full width half maximum of a 1-d array
    Returns the indices of the array elements left and right closest to the
    maximum that cross below half the maximum value
    """

    assert (type(array) == type(np.ndarray([])))

    # Find the maximum in the interior of the array
    fwhm = 0.5 * array[1:-1].max()
    max_idx = array[1:-1].argmax() + 1
    # Divide the intervall in halfs at the peak and find the index of the
    # value in the left half of the intervall before it increases over max/2

    # The FWHM is between the indices closest to the maximum, whose values
    # are below 0.5 times the max
    try:
        left_idx = np.argwhere(array[1: max_idx] < fwhm)[-1] + 1
    except IndexError:
        # This can occurs if there is no value in the array smaller than
        # 0.5*max.
        # In this case, return the left limit of the array as the FWHM.
        left_idx = 0
    try:
        right_idx = np.argwhere(array[max_idx:-1] < fwhm)[0] + max_idx
    except IndexError:
        # Subtract one because the index of an array with size n is
        # within 0..n-1
        right_idx = array.size() - 1

    return np.array([int(left_idx), int(right_idx)])


class gauss_fixed_mu(object):
    def __init__(self, mu):
        self.mu = mu
    def __call__(self, x, sigma):
        return np.exp(-0.5 * (x - self.mu) * (x - self.mu) / (sigma * sigma))


def width_gaussian(y, mu, imin=4, imax=8, l2min=4):
    """
    Fit a Gaussian on y, centered around mu
    on the interval [mu - i:mu - i], where i = imin..imax

    The best fit is given by the one, which minimizes the L2 norm on the
    intervall [mu - l2min:mu + l2min]

    Return the width of the best fit

    Input:
            y  : Data to fit
            mu : Center of the gauss

    Output:
            sigma:  Width of the resulting fit
            i    :  i value that achieved the best fit
            err  :  L2 error of the fit
    """

    assert(l2min <= imin)

    gaussian_fun = gauss_fixed_mu(mu)
    def err_fun(p, y, x):
        return (y - gaussian_fun(x, p))

    irange = np.arange(imin, imax, dtype='int')
    xrg = np.arange(0, y.shape[0])
    p = np.zeros(irange.shape[0])
    L2 = np.zeros(irange.shape[0])

#    plt.figure()
#    plt.plot([mu], [1.0], 'ko')

    for idx, i in enumerate(irange):
        idx_lo = int(max(mu - i, 0))
        idx_up = int(min(mu + i, y.shape[0]))
        fit_range = xrg[idx_lo:idx_up]
        fit_data = y[idx_lo:idx_up] / y[idx_lo:idx_up].max()

        p0 = [1.0]
        p_i, success_i = leastsq(err_fun, p0, args=(fit_data, fit_range), maxfev=1000)
        p[idx] = p_i

        mid_idx = fit_data.shape[0] / 2
        l2_idx = np.arange(mid_idx - l2min, mid_idx + l2min, 1, dtype='int')
        L2[idx] = ((fit_data[l2_idx] - gaussian_fun(p_i, fit_range[l2_idx])) ** 2.0).sum() / float(l2_idx.shape[0])

#        print i, p[idx], L2[idx], success_i
#
#        label_str = 'sigma = %5.2f, L2 = %5.2f' % (p[idx], L2[idx])
#        plt.plot(fit_range, fit_data)
#        plt.plot(fit_range, gaussian_fun(fit_range, p[idx]), label=label_str)

    best_fit_idx = L2.argmin()
#    idx_lo = int(max(mu - irange[best_fit_idx], 0))
#    idx_up = int(min(mu + irange[best_fit_idx], y.shape[0]))
#
#    plt.figure()
#    plt.plot([mu], [1.0], 'ko')
#    fit_range = xrg[idx_lo:idx_up]
#    fit_data = y[idx_lo:idx_up] / y[idx_lo:idx_up].max()
#    plt.plot(fit_range, fit_data)
#    plt.plot(fit_range, gaussian_fun(fit_range, p[best_fit_idx]), 'k', lw=3)

#    plt.legend()
#    plt.show()
    return (p[best_fit_idx], irange[best_fit_idx], L2[best_fit_idx])


class TrackingError(Exception):
    pass


def find_closest_region(frame, thresh_amp, x0, max_dist=2.0, verbose=False):
    """
    Returns the contiguous area above a threshold in a frame, whose centroid coordinates 
    are closest to x0

    Input:
        frame       : Input frame
        thresh_amp  : Threshold for area separation
        x0          : Distance to centroid
        max_dist    : Maximal distance to centroid

    Output:
        Binary image with contiguous area marked with 1
    """

    # Label all contiguous regions with ore than 60% of the original intensity
    labels = pm.label(frame > thresh_amp)
    # Get the area of all contiguous regions
    blob_area = pm.blob(labels, 'area', output='data')
    # Get the controid of all contiguous regions
    blob_cent = pm.blob(labels, 'centroid', output='data')
    if (verbose):
        print 'x0 = (%f, %f)' % (x0[0], x0[1])
        print 'Labelling found %d regions: ' % labels.max()
        for i in np.arange(labels.max()):
            print 'Region: %d, centroid at %d, %d, area: %d' % (i, blob_cent[i, 1], blob_cent[i, 0], blob_area[i])

    if (blob_cent.size < 1):
        raise TrackingError

    # We now have a bunch of contiguous regions.
    # Loop over the regions that are at least 10% of the largest region
    # and find the one, whose centroid is closest to the  last known position 
    # of the blob

    min_idx = -1   
    min_dist_frame = np.sqrt(frame.shape[0] * frame.shape[1])     # Maximal distance on a 64x64 grid

    for d_idx, i in enumerate(blob_area):
        # Compute distance of current areas centroid to the last centroids position
        dist = np.sqrt((blob_cent[d_idx, 1] - x0[1]) ** 2 +
                       (blob_cent[d_idx, 0] - x0[0]) ** 2)
        if (verbose):
            print 'Region %d, center: x=%d, y=%d, A=%f, distance to last centroid: %f' %\
                (d_idx, blob_cent[d_idx, 0], blob_cent[d_idx, 1], i, dist)

        # Skip areas who are less than 10% of the original
        if (i < 0.1 * blob_area.max()):
            if(verbose):
                print 'passing blob with area %f, d_idx = %d' % (i, d_idx)
            continue

        if (dist < min(max_dist, min_dist_frame)):
            min_dist_frame = dist
            min_idx = d_idx
            if (verbose):
                print 'Accepted'

    # If min_dist_frame is still sqrt(8192), no area was selected and the
    # blob could not be tracked successfully
    if (min_idx is -1):
        print 'No peak satisfying criteria.'
        raise TrackingError

    x_centroid = blob_cent[min_idx]

    # Compute the x and y COM coordinates of the blob, store
    blob_mask = labels != (min_idx + 1)
    event_masked = np.ma.MaskedArray(frame,
                                     mask=blob_mask,
                                     fill_value=0)
    #plt.figure()
    #plt.subplot(131)
    #plt.contourf(labels)
    #plt.colorbar()

    #plt.subplot(132)
    #plt.contourf(frame, 64)
    #plt.colorbar()

    #plt.subplot(133)
    #plt.contourf(event_masked)
    #plt.colorbar()

    #plt.show()

    return (x_centroid, event_masked)


def tracker(frames, x0, event, thresh_amp, thresh_dist, blob_ext,
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
        amp[tau] = event_frame[z0_last, R0_last]

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

    if (plots):
        plt.show()

    return tau, amp, xycom, xymax, width_rad, width_pol


def find_sol_pixels(shotnr, frame_info=None, rz_array=None,
                    datadir='/Users/ralph/source/blob_tracking/test_data'):
    """
    Returns the indices of the pixels in between the separatrix and the LCFS.
    """

    s = readsav('%s/separatrix.sav' % (datadir), verbose=False)

    gap_idx_mask = ((s['rmid'].reshape(64, 64) > s['rmid_sepx']) &
                    (s['rmid'].reshape(64, 64) < s['rmid_lim']))

    return np.argwhere(gap_idx_mask)


def find_sol_mask(shotnr, frame_info=None, rz_array=None,
                  datadir='/Users/ralph/source/blob_tracking/test_data'):
    """
    Returns a mask for the pixels in between the separatrix and the LCFS.
    """
    s = readsav('%s/separatrix.sav' % (datadir), verbose=False)

    return ((s['rmid'].reshape(64, 64) > s['rmid_sepx']) &
            (s['rmid'].reshape(64, 64) < s['rmid_lim']))


def blob_in_sol(trail, good_domain, logger=None):
    """
    Returns a bool array of the indices, in which the COM of a blob is
    in the SOL (good_domain)
    """

    try:
        # Consider only the positions, where the blob is in the good domain
        blob_pos = trail.get_trail_com()
        good_pos_idx = np.array([i in good_domain for i in
                                 blob_pos.round().astype('int').tolist()])

    except:
        good_pos_idx = np.ones_like(trail.get_tau())
        if (logger is not None):
            logger.info('This should not happen. Ignoring trigger domain')

    good_pos_idx = good_pos_idx[:-1]
    return good_pos_idx


if __name__ == "__main__":
    print "helper_functions.py"

# End of file helper_functions.py