#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib as mpl
from figure_defs import set_rcparams_pres_small
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt

"""
Plot mean COM velocities, amplitudes, and gasss at maximum, weighted mean. 
Values gathered for blob evolution in SOL only.

Using a CA threshold of 2.0

"""

save = False
filename = '/Users/ralph/source/blob_tracking/scan_120217/results_sol_thresh_15.eps'
thresh = 2.0

shotlist = np.array([1100811007, 1100811009, 1100811010, 1100811011, 1100811013])
neng     = np.array([0.16, 0.21, 0.21, 0.28, 0.33])

nblobs     = np.array([189, 94, 175, 180, 132])

vrad       = np.array([188.80, 248.54, 211.57, 229.81, 247.70])
err_vrad   = np.array([103.86, 141.39, 136.51, 151.46, 159.14])

vpol       = np.array([-29.39, -102.23, -125.01, -54.39, -87.58])
err_vpol   = np.array([188.34, 216.01, 239.51, 221.36, 300.32])

ell_rad    = np.array([0.57, 0.59, 0.59, 0.58, 0.57])
err_ellrad = np.array([0.04, 0.06, 0.06, 0.06, 0.07])
ell_pol    = np.array([0.55, 0.55, 0.57, 0.56, 0.54])
err_ellpol = np.array([0.05, 0.04, 0.05, 0.05, 0.05])

corr = np.array([-0.12, 0.08, -0.11, -0.19, -0.11])


np.savez('scan_100811/results_ds_sol_gauss_thresh20.npz', thresh = thresh, shotlist = shotlist, \
neng = neng, nblobs = nblobs, vrad = vrad, err_vrad = err_vrad, vpol = vpol, err_vpol = err_vpol,\
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
