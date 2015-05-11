#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib as mpl
#mpl.use('AGG')
import matplotlib.pyplot as plt
from misc.phantom_helper import make_rz_array
#from load_mdsframes import load_mdsframes
from scipy.interpolate import griddata
from scipy.stats import kurtosis, skew
#import idlsave
from scipy.io.idl import readsav

"""
Compute statistics for single pixels of the GPI array
"""


#shots = [ 1100803009, 1100803011, 1100803012, 1100803012,\
#    1100803015, 1100803020]
frame0 = 0
shots = [1111208020]

datadir = '/Users/ralph/uni/cmod/tmp_data/phantom/data/'
datadir2 = '/Users/ralph/source/blob_tracking/'
save_plots = False

for shotnr in shots:
#    frames, frame_info = load_mdsframes( shotnr )
    datafile = np.load('%s/%d_frames_normalized.npz' % (datadir, shotnr),
                       mmap_mode='r')
    print datafile.keys()
    frames = datafile['frames_normalized_rms']
    frame_info = datafile['frame_info']
    frames_mean = datafile['frames_mean']
    frames_rms = datafile['frames_rms']
#    frames = frames[frame0:, :, :]

    print np.shape(frames)
    try:
        s = readsav('%d/%d_separatrix.sav' % (shotnr, shotnr),
                    verbose=False)
    except:
        print 'Cannot load separatrix data for this shot.'
        print 'Using default values for #1120217008'
        s = readsav('%s/1120217008/1120217008_separatrix.sav' % (datadir2),
                    verbose=False)

    rz_array, transform_data = make_rz_array(frame_info)
    print np.shape(rz_array)

    r_flat = rz_array[:, :, 0].reshape(64 * 64)
    z_flat = rz_array[:, :, 1].reshape(64 * 64)

    # define grid.
    xxi, yyi = np.meshgrid(np.linspace(rz_array[:, :, 0].min(),
                                       rz_array[:, :, 0].max(), 64),
                           np.linspace(rz_array[:, :, 1].min(),
                                       rz_array[:, :, 1].max(), 64))
    xyi = np.concatenate((xxi[:, :, np.newaxis], yyi[:, :, np.newaxis]),
                         axis=2)
    # grid the data.
    zi = griddata(rz_array.reshape(64 * 64, 2),
                  frames[666, :, :].reshape(64 * 64), xyi.reshape(64 * 64, 2),
                  method='linear')
    xxi_flat, yyi_flat = xxi.reshape(64 * 64), yyi.reshape(64 * 64)

    print 'shotnr %d' % shotnr
    print 'Separatrix at %f cm' % (s.rmid_sepx * 100.)
    print 'Limiter shadow at %f cm' % (s.rmid_lim * 100.)

    # Flat indices of points between separatrix and limiter shadow
    gap_idx = np.all(((xxi_flat > s.rmid_sepx * 100),
                      (xxi_flat < s.rmid_lim * 100)), axis=0)
    # 2d array indices of pixels between separatrix and limiter shadow
    gap_idx_unr = np.unravel_index(np.where(gap_idx)[0], (64, 64))
    gap_idx_unr = np.concatenate((gap_idx_unr[0][:, np.newaxis],
                                  gap_idx_unr[1][:, np.newaxis]),
                                 axis=1)

    # Analyze only pixel at poloidal index 32
    n = np.size(np.where(gap_idx_unr[:, 0] == 32)[0])
    # Take only 5 positions
    r_idx = np.where(gap_idx_unr[:, 0] == 32)[0][::np.round(n / 5) + 1]
    n_r = np.size(np.array(r_idx))

    # For frames, axis1 is the poloidal direction,
    # axis2 is the radial direction
    pixels = [frames[:, gap_idx_unr[r, 0], gap_idx_unr[r, 1]] for r in r_idx]
    # For rz_array, axis0 is the radial pixel,
    # axis1 is the poloidal pixel
    radii = [rz_array[gap_idx_unr[r, 0], gap_idx_unr[r, 1], 0] for r in r_idx]
    pol = [rz_array[gap_idx_unr[r, 0], gap_idx_unr[r, 1], 1] for r in r_idx]

    # Take fluctuation and RMS from pre-computed array
    rms = np.array([frames_rms[gap_idx_unr[r, 0],
                               gap_idx_unr[r, 1]] for r in r_idx])
    mean = np.array([frames_mean[gap_idx_unr[r, 0],
                                 gap_idx_unr[r, 1]] for r in r_idx])

    # Store fluctuations/mean, i_mean, skew and flatness for
    # each radial position
    skew_arr = np.zeros(n_r)
    kurt_arr = np.zeros(n_r)

    # Interpolate mean and rms on the new grid
    frames_mean_ip = griddata(rz_array.reshape(64 * 64, 2),
                              frames_mean.reshape(64 * 64),
                              xyi.reshape(64 * 64, 2),
                              method='linear').reshape(64, 64)
    frames_rms_ip = griddata(rz_array.reshape(64 * 64, 2),
                             frames_rms.reshape(64 * 64),
                             xyi.reshape(64 * 64, 2),
                             method='linear').reshape(64, 64)

    print np.shape(frames_mean_ip)
    print np.shape(frames_rms_ip)
    print np.shape(frames_rms_ip / frames_mean_ip)

    f = plt.figure()
    plt.title('#%d, $\\bar{I}$' % shotnr)
    plt.contour(xyi[:, :, 0], xyi[:, :, 1], frames_mean_ip, 15,
                linewidth=0.5, colors='k')
    plt.contourf(xyi[:, :, 0], xyi[:, :, 1], frames_mean_ip, 15,
                 cmap=plt.cm.hot)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')
    plt.colorbar()
    if save_plots:
        f.savefig('%d/statistics/%d_mean.eps' % (shotnr, shotnr))

    f = plt.figure()
    plt.title('#%d, $I_{RMS}$' % shotnr)
    plt.contour(xyi[:, :, 0], xyi[:, :, 1], frames_rms_ip,
                15, linewidth=0.5, colors='k')
    plt.contourf(xyi[:, :, 0], xyi[:, :, 1], frames_rms_ip,
                 15, cmap=plt.cm.hot)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')
    plt.colorbar()
    if save_plots:
        f.savefig('%d/statistics/%d_rms.eps' % (shotnr, shotnr))

    f = plt.figure()
    plt.title('#%d, RMS over mean' % shotnr)
    plt.contour(xyi[:, :, 0], xyi[:, :, 1], (frames_rms_ip / frames_mean_ip),
                15, linewidth=0.5, colors='k')
    plt.contourf(xyi[:, :, 0], xyi[:, :, 1], (frames_rms_ip / frames_mean_ip),
                 15, cmap=plt.cm.hot)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')
    plt.colorbar()
    if save_plots:
        f.savefig('%d/statistics/%d_rms_over_mean.eps' % (shotnr, shotnr))

    plt.figure()
    for idx, pix in enumerate(pixels):
        print 'Statistics at R,z = %4.3f cm, %3.2f cm' % (radii[idx], pol[idx])
        plt.subplot(2, 1, 1)
        plt.plot(pix[::20] + 4 * idx, label='R=%4.2f cm' % (radii[idx]))
        plt.legend(loc='upper left')

        plt.subplot(2, len(pixels), len(pixels) + idx + 1)
    #    n, bins, patches = plt.hist(pix, bins=25, log=True,
    #                                normed=True, histtype='step', lw=2)
        n, bins, patches = plt.hist(pix, bins=25, log=True)

        skew_arr[idx] = skew(pix)
        kurt_arr[idx] = kurtosis(pix)

    if save_plots:
        F = plt.gcf()
        F.savefig('%d/statistics/%d_pixel_timeseries.png' % (shotnr, shotnr),
                  dpi=300)
    plt.close()

    F = plt.figure()
    F.text(0.5, 0.95, 'shot #%d' % shotnr, ha='center')
    plt.subplot(221)
    plt.plot(radii, mean, 'ko-')
    plt.xlabel('R / cm')
    plt.ylabel('$\\bar{I}$')

    plt.subplot(222)
    plt.plot(radii, rms / mean, 'ko-')
    plt.xlabel('R / cm')
    plt.ylabel('$I_{RMS} / \\bar{I}$')

    plt.subplot(223)
    plt.plot(radii, skew_arr, 'ko-')
    plt.xlabel('R / cm')
    plt.ylabel('Skewness')

    plt.subplot(224)
    plt.plot(radii, kurt_arr, 'ko-')
    plt.xlabel('R / cm')
    plt.ylabel('Kurtosis')

    if save_plots:
        F.savefig('%d/statistics/%d_pixel_statistics.eps' % (shotnr, shotnr))
        plt.close()

    plt.show()
