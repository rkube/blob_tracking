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

Using a CA threshold of 1.5
Gaussian width detection, use weighted mean for each blob

"""

save = False
filename = '/Users/ralph/source/blob_tracking/scan_120217/results_sol_thresh_15_gauss_wmean.eps'
thresh = 1.5

#shotlist = [1120217008, 1120217010, 1120217011, 1120217012, 1120217014, 1120217015, 1120217016, \
#    1120217019, 1120217020, 1120217021]

shotlist = [1120217008, 1120217010, 1120217011, 1120217015, 1120217016, 1120217017]


#neng     = np.array([0.13, 0.1395, 0.1405, 0.32, 0.33, 0.34, 0.40, 0.48])
neng       = np.array([0.133, 0.14, 0.14, 0.19, 0.21, 0.29, 0.30, 0.34, 0.40, 0.48]) 

nblobs     = np.array([499, 397, 436, 4, 32, 3, 4, 193, 274, 295])

vrad       = np.array([183.22, 154.18, 161.29, 123.50, 294.45, 294.71, 287.43, 355.53, 415.53, 411.50])
err_vrad   = np.array([141.03, 155.02, 143.04, 104.44, 240.16, 86.40, 73.22, 198.17, 242.69, 251.38])

vpol       = np.array([-152.42, -149.45, -116.60, -160.95, -155.01, -369.05, -195.17, -104.57, -83.56, -81.90])
err_vpol   = np.array([263.13, 324.34, 311.13, 338.92, 236.98, 169.04, 263.79, 241.45, 313.64, 346.85])

ell_rad    = np.array([0.67, 0.68, 0.67, 0.77, 0.71, 0.72, 0.66, 0.71, 0.72, 0.00])
err_ellrad = np.array([0.08, 0.07, 0.08, 0.04, 0.08, 0.14, 0.10, 0.13, 0.16, 0.00])
ell_pol    = np.array([0.65, 0.66, 0.65, 0.73, 0.68, 0.69, 0.66, 0.68, 0.70, 0.69])
err_ellpol = np.array([0.08, 0.07, 0.07,0.05, 0.07, 0.11, 0.10, 0.09, 0.08, 0.09])

corr = np.array([-0.09, -0.14, -0.20, -0.12, -0.10, 0.55, -0.86, -0.28, -0.22, -0.11])

np.savez('scan_120217/results_ds_sol_gauss_thresh15_0418.npz', thresh = thresh,  shotlist = shotlist,\
neng = neng, nblobs = nblobs, vrad = vrad, err_vrad = err_vrad, vpol = vpol, err_vpol = err_vpol, \
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
