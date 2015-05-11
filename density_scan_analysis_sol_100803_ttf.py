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
Only blobs in SOL
Gaussian width detection

"""

save = False
filename = '/Users/ralph/source/blob_tracking/scan_100803/results_sol_thresh_25_ttf.eps'


shotlist = np.array([1100803005, 1100803006, 1100803007, 1100803008, 1100803009, 1100803011, 1100803012, 1100803013, 1100803015, 1100803020 ] )
neng     = np.array([0.151, 0.14, 0.149, 0.20, 0.18, 0.239, 0.241, 0.28, 0.309, 0.311] )

nblobs     = np.array([135, 118, 155, 150, 147, 118, 110, 132, 155, 155])

vrad       = np.array([350.24, 339.61, 303.22, 310.04, 304.69, 377.13, 382.88, 407.54, 426.86, 391.42])
err_vrad   = np.array([142.10, 153.27, 160.28, 154.79, 168.03, 168.33, 191.67, 199.76, 216.46, 202.23])

vpol       = np.array([-141.71, -142.25, -117.67, -76.20, -89.83, -107.08, -91.17, -80.67, -85.29,  -50.39])
err_vpol   = np.array([228.60, 243.70, 251.28, 241.04, 226.45, 243.83, 259.00, 297.00, 257.64, 293.16])

ell_rad    = np.array([0.52, 0.52, 0.52, 0.52, 0.52, 0.52, 0.55, 0.55, 0.54, 0.53])
err_ellrad = np.array([0.04, 0.04, 0.06, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05]) 
ell_pol    = np.array([0.55, 0.56, 0.57, 0.59, 0.58, 0.59, 0.61, 0.58, 0.56, 0.58])
err_ellpol = np.array([0.08, 0.08, 0.07, 0.09, 0.09, 0.08, 0.06, 0.07, 0.06, 0.09])

corr = np.array([-0.44, -0.46, -0.24, -0.24, -0.36, -0.13, -0.34, -0.26, -0.28, -0.22])

np.savez('scan_100803/results_ds_sol_gauss_thresh25_ttf.npz', thresh = 2.5, shotlist = shotlist, \
neng = neng, nblobs = nblobs, vad = vrad, err_vrad = err_vrad, vpol = vpol, err_vpol = err_vpol, \
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
