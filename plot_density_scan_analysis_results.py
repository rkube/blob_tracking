#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib as mpl
from plotting.figure_defs import set_rcparams_pres_small
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt



datafile_names = ['scan_100803/results_ds_sol_gauss_thresh25_0418.npz', 'scan_120217/results_ds_sol_gauss_thresh15_0418.npz',\
'scan_100811/results_ds_sol_gauss_thresh20.npz']
symbols = ['bv', 'gs', 'r^']


fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize = (12,6.5))
#fig.text(.5,.95, 'Shots 1100803005-020', ha='center')
ax = axs[0,0]
ax.set_ylabel('Radial velocity / m/s', fontsize = 24)
ax.set_ylim( (-100, 800) )

ax = axs[0,1]
ax.set_ylabel('Poloidal velocity / m/s', fontsize = 24)


ax = axs[1,0]
ax.set_xlabel('$n_e / n_G$', fontsize=24)
ax.set_ylabel('Radial  length / cm', fontsize = 24)
ax.set_ylim( (0.4, 1.0))

#plt.subplot(224)
ax = axs[1,1]
ax.set_xlabel('$n_e / n_G$', fontsize=24)
ax.set_ylabel('Poloidal  length / cm', fontsize = 24)

ax.set_xticks([0.15, 0.2, 0.25, 0.3, 0.4, 0.5])
ax.set_xlim((0.13, 0.50))

plot_symbols = ['o', 's', 'v']
color_idx = [ 'b', 'g', 'r' ]
legends = ['B=4.0T', 'B=5.4T', 'B=8.0T']


for idx, fn in enumerate(datafile_names):
    datafile = np.load(fn)
    thresh = datafile['thresh']
    shotlist = datafile['shotlist']
    neng = datafile['neng']
    nblobs = datafile['nblobs']
    vrad = datafile['vrad']
    err_vrad = datafile['err_vrad']
    vpol = datafile['vpol']
    err_vpol = datafile['err_vpol']
    ell_rad = datafile['ell_rad']
    err_ellrad = datafile['err_ellrad']
    ell_pol = datafile['ell_pol']
    err_ellpol = datafile['err_ellpol']
    corr= datafile['corr']


    ax = axs[0,0]
    ax.errorbar( neng, vrad, yerr = err_vrad, ls = ' ', color = color_idx[idx], marker = plot_symbols[idx], markersize = 8, label = legends[idx])
    
    ax = axs[0,1]
    ax.errorbar( neng, vpol, yerr = err_vpol, ls = ' ', color = color_idx[idx], marker = plot_symbols[idx], markersize = 8, label = legends[idx])
    
    ax = axs[1,0]
    ax.errorbar( neng, ell_rad, yerr = err_ellrad, ls = ' ', color = color_idx[idx], marker = plot_symbols[idx], markersize = 8, label = legends[idx])

    ax = axs[1,1]
    ax.errorbar( neng, ell_pol, yerr = err_ellpol, ls = ' ', color = color_idx[idx], marker = plot_symbols[idx], markersize = 8, label = legends[idx])
    
axs[0,0].legend(loc = 'lower right')
axs[0,1].legend(loc = 'lower right')
axs[1,0].legend(loc = 'lower right')
axs[1,1].legend(loc = 'lower right')

    
plt.show()
