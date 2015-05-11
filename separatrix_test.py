#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from misc.phantom_helper import make_rz_array
from scipy.interpolate import griddata
from plotting.separatrix_line import surface_line
#import idlsave
#from scipy.io.idl import readsav
from scipy.io import readsav

shotnr = 1111208020
datadir = '/Volumes/My Book Thunderbolt Duo/cmod_data/phantom/data/'
#datadir = '.'
df_frames = '%s/%10d_testframes.npz' % (datadir, shotnr)
#df_frames = '%s/%10d/%10d_testframes.npz' % (datadir, shotnr, shotnr)
df_separatrix = '%s/%10d_separatrix.sav' % (datadir, shotnr)
#df_separatrix = '%s/%10d/%10d_separatrix.sav' % (datadir, shotnr, shotnr)

datafile = np.load(df_frames)
frames = datafile['frames']
#frames = datafile['frames_normalized_rms']
frame_info = datafile['frame_info']
#s = readsav('1100803005/1100803005_separatrix.sav', verbose=False)
s = readsav(df_separatrix, verbose=True)

# Find all flux surfaces in the gap between separatrix and limiter shadow in
# the original, non-interpolated image
gap_idx_mask = np.invert((s['rmid'].reshape(64, 64) > s['rmid_sepx']) &
                         (s['rmid'].reshape(64, 64) < s['rmid_lim']))
gap_idx_mask_inv = (s['rmid'].reshape(64, 64) > s['rmid_sepx']) &\
                   (s['rmid'].reshape(64, 64) < s['rmid_lim'])
frame_m = ma.array(frames[666, :, :], mask=gap_idx_mask)


separatrix_pxs = surface_line(s['rmid'].reshape(64, 64) >
                              s['rmid_sepx'], mode='max')
limiter_pxs = surface_line(s['rmid'].reshape(64, 64) <
                           s['rmid_lim'], mode='min')

rz_array, transform_data = make_rz_array(frame_info)

r_flat = rz_array[:, :, 0].reshape(64 * 64)
z_flat = rz_array[:, :, 1].reshape(64 * 64)

# define grid.
xxi, yyi = np.meshgrid(np.linspace(rz_array[:, :, 0].min(),
                                   rz_array[:, :, 0].max(), 64),
                       np.linspace(rz_array[:, :, 1].min(),
                                   rz_array[:, :, 1].max(), 64))
xyi = np.concatenate((xxi[:, :, np.newaxis],
                      yyi[:, :, np.newaxis]), axis=2)
print xyi.shape, xyi[0, 0, :], xyi[63, 63, :]
print rz_array.shape, rz_array[0, 0, :], rz_array[63, 63, :]
# grid the data.
zi = griddata(rz_array.reshape(64 * 64, 2), frames[666, :, :].reshape(64 * 64),
              xyi.reshape(64 * 64, 2), method='linear')

plt.figure()
plt.title('#%10d: phantom view' % (shotnr))
plt.contourf(zi.reshape(64, 64), 32, cmap=plt.cm.hot)
plt.colorbar()
xxi_flat = xxi.reshape(64 * 64)
yyi_flat = yyi.reshape(64 * 64)

fig = plt.figure()
plt.title('#%10d: phantom view' % (shotnr))
plt.contourf(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64), 32,
             cmap=plt.cm.hot)
plt.colorbar()
plt.xlabel('R / cm')
plt.ylabel('Z / cm')

plt.plot((xyi[:, :, 0])[gap_idx_mask], (xyi[:, :, 1])[gap_idx_mask],
         'ko', markersize=1.5)

print 'Trigger box: lower left: R=%f Z=%f, width = %f, height = %f' % \
    (xyi[16, 40, 0], xyi[16, 40, 1], (xyi[16, 50, 0] - xyi[16, 40, 0]),
    (xyi[48, 40, 1] - xyi[16, 40, 1]))
tb_lower_left = (xyi[16, 40, 0], xyi[16, 40, 1])
tb_width = (xyi[16, 50, 0] - xyi[16, 40, 0])
tb_height = (xyi[48, 40, 1] - xyi[16, 40, 1])

# Plot the triggering domain
triggering_box = mpatches.Rectangle(tb_lower_left, width=tb_width,
                                    height=tb_height, fill=False,
                                    ls='dashdot', ec='w', lw=3)
ax = fig.gca()
ax.add_patch(triggering_box)

# Plot the separatrix
sep_x = [xyi[i, separatrix_pxs[i], 0] for i in np.arange(64)]
sep_y = [xyi[i, separatrix_pxs[i], 1] for i in np.arange(64)]
plt.plot(sep_x, sep_y, 'w--', linewidth=4)

lim_x = [xyi[i, limiter_pxs[i], 0] for i in np.arange(64)]
lim_y = [xyi[i, limiter_pxs[i], 1] for i in np.arange(64)]
plt.plot(lim_x, lim_y, 'w-.', linewidth=4)

print 'Separatrix at %f cm' % (s.rmid_sepx * 100.)
print 'Limiter shadow at %f cm' % (s.rmid_lim * 100.)

plt.show()
