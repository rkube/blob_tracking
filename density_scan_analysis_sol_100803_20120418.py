#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib as mpl
from figure_defs import set_rcparams_pres_small
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt

"""
Plot mean COM velocities, amplitudes, and FWHM at maximum. Values gathered for blob evolution
in SOL only.

Using a CA threshold of 2.0
Only blobs in SOL
Gaussian width detection, use weighted mean for each blob

"""

save = False
filename = '/Users/ralph/source/blob_tracking/scan_100803/results_sol_thresh_25_ttf.eps'


shotlist = np.array([1100803005, 1100803006, 1100803007, 1100803008, 1100803009, 1100803011, 1100803012, 1100803013, 1100803015, 1100803020 ] )
neng     = np.array([0.151, 0.14, 0.149, 0.20, 0.18, 0.239, 0.241, 0.28, 0.309, 0.311] )

nblobs     = np.array([499, 517, 567, 590, 570, 512, 501, 502, 509, 532])

vrad       = np.array([284.66, 271.60, 267.72, 269.54, 274.83, 318.07, 334.00, 384.37, 395.05, 384.51])
err_vrad   = np.array([152.06, 159.87, 151.53, 154.26, 163.28, 164.02, 170.32, 199.29, 219.84, 215.06])

vpol       = np.array([-134.24, -111.39, -112.90, -55.52, -78.13, -52.03, -63.99, -53.87, -53.02, -36.80])
err_vpol   = np.array([264.83, 245.22, 270.93, 249.34, 261.52, 257.92, 257.02, 294.12, 286.64, 297.77])

ell_rad    = np.array([0.56, 0.56, 0.56, 0.56, 0.55, 0.57, 0.58, 0.58, 0.58, 0.57])
err_ellrad = np.array([0.06, 0.06, 0.07, 0.06, 0.07, 0.06, 0.06, 0.06, 0.06, 0.06])
ell_pol    = np.array([0.55, 0.56, 0.55, 0.55, 0.54, 0.56, 0.56, 0.56, 0.56, 0.55])
err_ellpol = np.array([0.06, 0.06, 0.06, 0.06, 0.07, 0.05, 0.05, 0.06, 0.05, 0.06])

corr = np.array([-0.36, -0.47, -0.27, -0.27, -0.29, -0.28, -0.25, -0.26, -0.20, -0.22])


np.savez('scan_100803/results_ds_sol_gauss_thresh25_0418.npz', thresh = 2.5, shotlist = shotlist, \
neng = neng, nblobs = nblobs, vrad = vrad, err_vrad = err_vrad, vpol = vpol, err_vpol = err_vpol, \
ell_rad = ell_rad, err_ellrad = err_ellrad, ell_pol = ell_pol, err_ellpol = err_ellpol, corr = corr)

#np.savez('scan_100803/results_blob_tracking_gauss.npz', shotlist = shotlist, neng = neng, 
#nblobs = nblobs, vrad = vrad, err_vard = err_vrad, vpol = vpol, err_vpol = err_vpol,
#ell_rad = ell_rad, err_ellrad = err_ellrad, ell_pol = ell_pol, err_ellpol = err_ellpol,
#corr = corr)

f = plt.figure()
ax1 = plt.subplot(211)
plt.plot( neng, nblobs, 'ko', )
plt.ylabel('# of blob events')

ax2 = plt.subplot(212, sharex = ax1)
plt.plot( neng, corr, 'gs' )
plt.ylabel('Correlation between v_r and $\ell_{pol}$')
plt.xlim( ( 0.98*neng.min(), 1.02*neng.max() ) )

plt.xlabel('$n_e / n_G$')

fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize = (12,6.5))
#fig.text(.5,.95, 'Shots 1100803005-020', ha='center')
ax = axs[0,0]
ax.errorbar( neng, vrad, yerr = err_vrad, ls = ' ', color = 'g', marker = '>', markersize = 8 )
ax.set_ylabel('Radial velocity / m/s', fontsize = 16)

ax = axs[0,1]
ax.errorbar( neng, vpol, yerr = err_vpol, ls = ' ', color = 'k', marker = 'v', markersize = 8 )
ax.set_ylabel('Poloidal velocity / m/s', fontsize = 16)


ax = axs[1,0]
ax.errorbar( neng, ell_rad, yerr = err_ellrad, ls = ' ', color = 'g', marker = '>', markersize = 8 )
ax.set_xlabel('$n_e / n_G$', fontsize=24)
ax.set_ylabel('Radial  length / cm', fontsize = 16)

#plt.subplot(224)
ax = axs[1,1]
ax.errorbar( neng, ell_pol, yerr = err_ellpol, ls = ' ', color = 'k', marker = 'v', markersize = 8 )
ax.set_xlabel('$n_e / n_G$', fontsize=24)
ax.set_ylabel('Poloidal  length / cm', fontsize = 16)

ax.set_xticks([0.15, 0.2, 0.25, 0.3])
ax.set_xlim((0.13, 0.32))

if ( save == True ):
    f.savefig(filename)

plt.show()
