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

Using a CA threshold of 2.5

Skipped shots: 1100803011 ( Error in blob tracking, peak 863 )

"""

save = False
filename = '/Users/ralph/source/blob_tracking/scan_100803/results_sol_thresh_25.eps'

shotlist = np.array([1100803005, 1100803006, 1100803007, 1100803008, 1100803009, 1100803011, 1100803012, 1100803013, 1100803015, 1100803020 ] )
neng     = np.array([0.151, 0.14, 0.149, 0.20, 0.18, 0.239, 0.241, 0.28, 0.309, 0.311] )

nblobs     = np.array([407, 401, 446, 451, 450, 445, 459, 487, 466, 437])

vrad       = np.array([305.81, 308.35, 303.59, 307.64, 315.72, 371.68, 390.93, 447.73, 460.97, 464.86])
err_vrad   = np.array([159.60, 158.27, 150.57, 159.72, 165.91, 178.84, 188.29, 208.61, 215.67, 227.22])

vpol       = np.array([-145.23, -135.80, -107.97, -72.86, -98.72, -79.01, -79.02, -33.21, -23.56, -31.06])
err_vpol   = np.array([246.85, 253.26, 237.32, 224.50, 242.15, 263.34, 279.87, 313.97, 287.65, 317.99])

ell_rad    = np.array([0.63, 0.63, 0.62, 0.62, 0.61, 0.63, 0.65, 0.65, 0.64, 0.63])
err_ellrad = np.array([0.11, 0.10, 0.12, 0.10, 0.11, 0.12, 0.12, 0.11, 0.11, 0.11])
ell_pol    = np.array([0.65, 0.65, 0.63, 0.65, 0.64, 0.66, 0.67, 0.64, 0.61, 0.63])
err_ellpol = np.array([0.12, 0.12, 0.11, 0.11, 0.12, 0.13, 0.12, 0.12, 0.10, 0.12])

corr = np.array([-0.62, -0.58, -0.44, -0.41, -0.47, -0.34, -0.29, -0.29, -0.24, -0.25])

np.savez('scan_100803/results_ds_sol_gauss_thresh25.npz', thresh = 2.5, shotlist = shotlist, \
neng = neng, nblobs = nblobs, vrad = vrad, err_vrad = err_vrad, vpol = vpol, err_vpol = err_vpol, \
ell_rad = ell_rad, err_ellrad = err_ellrad, ell_pol = ell_pol, err_ellpol = err_ellpol, corr = corr)


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
