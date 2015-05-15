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
from scipy.optimize import leastsq
from scipy.io import readsav



def frac_to_idx(frac, f_min, f_max, nbins):
    """
    Given a float within in the interval [f_min,f_max], return the index in
    which of nbins equidistant bins it belongs
    """
    return np.floor((frac - f_min) * float(nbins) /
                    (f_max - f_min)).astype('int')



def fwhm(array):
    """
    Computes the full width half maximum of a 1-d array
    Returns the indices of the array elements left and right closest to the
    maximum that cross below half the maximum value

    Input:
        array:    ndarray

    Output:
        idx_arr:  ndarray, Indices of the array left and right closest to the
                           maximum where array crosses below half its maximum
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

    idx_arr = np.array([int(left_idx), int(right_idx)])

    return idx_arr


class gauss_fixed_mu(object):
    """
    Gaussian function functor with fixed mean
    """
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

    # the index range on which we fit
    irange = np.arange(imin, imax, dtype='int')
    # corresponding x range, in pixel
    xrg = np.arange(0, y.shape[0])
    # Best fit for given index range
    p = np.zeros(irange.shape[0])
    # L2 norm of best fit on given index range
    L2 = np.zeros(irange.shape[0])

    #plt.figure()
    #plt.plot([mu], [1.0], 'ko')

    for idx, i in enumerate(irange):
        # Build interval to fit on, [mu - i : mu + i], stop at interval boundaries though
        idx_lo = int(max(mu - i, 0))
        idx_up = int(min(mu + i, y.shape[0]))
        # x range, in pixel
        fit_range = xrg[idx_lo:idx_up]
        # Data to fit on, normalized to unity
        fit_data = y[idx_lo:idx_up] / y[idx_lo:idx_up].max()
        #print 'mu = %d, i = %d, idx_lo = %d, idx_up = %d' % (mu, i, idx_lo, idx_up)
        #print 'fit_range = ', fit_range, ', fit_data  = ', fit_data

        # Fit Gaussian on interval
        p0 = [1.0]
        p_i, success_i = leastsq(err_fun, p0, args=(fit_data, fit_range), maxfev=1000)
        p[idx] = p_i

        # Compute L2 norm of the fit
        mid_idx = fit_data.shape[0] / 2
        l2_idx = np.arange(max(0, mid_idx - l2min), min(63, mid_idx + l2min), 1, dtype='int')

        #print 'L2.shape = ', L2.shape, 'idx = ', idx, ', fit_data.shape = ', fit_data.shape, ', l2_idx = ', l2_idx

        #L2[idx] = ((fit_data[l2_idx] - gaussian_fun(p_i, fit_range[l2_idx])) ** 2.0).sum() / float(l2_idx.shape[0])
        L2[idx] = ((fit_data - gaussian_fun(p_i, fit_range)) ** 2.0).sum() / float(fit_range.shape[0])

        label_str = 'sigma = %5.2f, L2 = %5.2f' % (p[idx], L2[idx])
        #plt.plot(fit_range, fit_data)
        #plt.plot(fit_range, gaussian_fun(fit_range, p[idx]), label=label_str)

    best_fit_idx = L2.argmin()
    idx_lo = int(max(mu - irange[best_fit_idx], 0))
    idx_up = int(min(mu + irange[best_fit_idx], y.shape[0]))

    #plt.figure()
    #plt.plot([mu], [1.0], 'ko')
    #fit_range = xrg[idx_lo:idx_up]
    #fit_data = y[idx_lo:idx_up] / y[idx_lo:idx_up].max()
    #plt.plot(fit_range, fit_data)
    #plt.plot(fit_range, gaussian_fun(fit_range, p[best_fit_idx]), 'k', lw=3)

    #plt.legend()
    #plt.show()
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



# End of file helper_functions.py
