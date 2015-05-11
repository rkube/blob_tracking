#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib as mpl
from plotting.figure_defs import set_rcparams_pres_small
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt

"""
Plot mean COM velocities, amplitudes, and FWHM at maximum. Values gathered for blob evolution
in SOL only.

Using a CA threshold of 2.5

"""

save = False
filename = '/Users/ralph/source/blob_tracking/scan_120217/results_sol_thresh_15.eps'


shotlist = np.array([120217008, 1120217010, 1120217011, 1120217017, 1120217018, \
    1120217019, 1120217020, 1120217021])
neng     = np.array([0.13, 0.1395, 0.1405, 0.32, 0.33, 0.34, 0.40, 0.48])

nblobs     = np.array([472, 384, 407, 90, 110, 91, 133, 170])

vrad       = np.array([177.29, 149.44, 158.15, 273.13, 348.81, 370.51, 429.73, 400.28])
err_vrad   = np.array([135.39, 148.58, 141.52, 212.74, 225.65, 193.61, 254.59, 253.04])

vpol       = np.array([-154.99, -162.51, -133.68, -106.80, -112.31, -138.84, -57.64, -84.35])
err_vpol   = np.array([266.99, 321.40, 305.07, 253.08, 241.40, 245.92, 297.73, 314.30])

ell_rad    = np.array([0.66, 0.67, 0.67, 0.69, 0.69, 0.71, 0.70, 0.70])
err_ellrad = np.array([0.09, 0.07, 0.08, 0.05, 0.08, 0.11, 0.07, 0.07])
ell_pol    = np.array([0.73, 0.76, 0.74, 0.71, 0.73, 0.70, 0.71, 0.70])
err_ellpol = np.array([0.12, 0.16, 0.11, 0.07, 0.17, 0.07, 0.07, 0.08])

corr = np.array([-0.19, -0.19, -0.23, -0.12, 0.01, -0.08, -0.26, -0.28])

np.savez('scan_120217/results_ds_sol_gauss_thresh15.npz', shotlist = shotlist, neng = neng, \
nblobs = nblobs, vrad = vrad, err_vard = err_vrad, vpol = vpol, err_vpol = err_vpol, \
ell_rad = ell_rad, err_ellrad = err_ellrad, ell_pol = ell_pol, err_ellpol = err_ellpol, corr = corr)


plt.figure()
ax1 = plt.subplot(211)
plt.plot( neng, nblobs, 'ko', )
plt.ylabel('# of blob events')

ax2 = plt.subplot(212, sharex = ax1)
plt.plot( neng, corr, 'gs' )
plt.ylabel('Correlation between v_r and $\ell_{pol}$')
plt.xlim( ( 0.98*neng.min(), 1.02*neng.max() ) )

plt.xlabel('$n_e / n_G$')

fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize = (12,6.5))
#fig.text(.5,.95, '#1120217008-021', ha='center')
ax = axs[0,0]
ax.errorbar( neng, vrad, yerr = err_vrad, ls = ' ', color = 'g', marker = '>', markersize = 8 )
ax.set_ylabel('Radial velocity / m/s', fontsize = 16)

ax = axs[0,1]
ax.errorbar( neng, vpol, yerr = err_vpol, ls = ' ', color = 'k', marker = 'v', markersize = 8 )
ax.set_ylabel('Poloidal velocity / m/s', fontsize = 16)

ax = axs[1,0]
ax.errorbar( neng, ell_rad, yerr = err_ellrad, ls = ' ', color = 'g', marker = '>', markersize = 8 )
ax.set_xlabel('$n_e / n_G$', fontsize = 24)
ax.set_ylabel('Radial  length / cm', fontsize = 16)

#plt.subplot(224)
ax = axs[1,1]
ax.errorbar( neng, ell_pol, yerr = err_ellpol, ls = ' ', color = 'k', marker = 'v', markersize = 8 )
ax.set_xlabel('$n_e / n_G$', fontsize = 24)
ax.set_ylabel('Poloidal  length / cm', fontsize = 16)

ax.set_xticks([0.15, 0.25, 0.35, 0.45])
ax.set_xlim((0.10, 0.50))

if ( save == True ):
    fig.savefig(filename)


plt.show()
