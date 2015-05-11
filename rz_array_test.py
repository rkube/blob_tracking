#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.projections import register_projection
#from phantom_helper import make_rz_array, GPI_projection
from misc.phantom_helper import make_rz_array
from scipy.interpolate import griddata
#import numpy.ma as ma
#from numpy.random import uniform, seed

shotnr = 1111208020
datadir = '/Users/ralph/uni/cmod/tmp_data/phantom/data'
f_file = '%s/%10d_frames_normalized.npz' % (datadir, shotnr)
datafile = np.load(f_file, mmap_mode='r')
frames = datafile['frames_normalized_rms']
frame_info = datafile['frame_info']

#register_projection(GPI_projection)

rz_array, transform_data = make_rz_array(frame_info)

print 'Shape of rz_array: ', rz_array.shape

#plt.figure()
#plt.title('R')
#plt.contourf(rz_array[:,:,0],32)
#plt.colorbar()
#
#plt.figure()
#plt.title('z')
#plt.contourf(rz_array[:,:,1],32)
#plt.colorbar()

plt.figure()
plt.title('test1')
plt.contour(frames[666, :, :], 15, linewidths=0.5, colors='k')
plt.contourf(frames[666, :, :], 15, cmap=plt.cm.hot)
plt.colorbar()

# define grid.
xxi, yyi = np.meshgrid(np.linspace(rz_array[:, :, 0].min(), rz_array[:, :, 0].max(), 64),
                       np.linspace(rz_array[:, :, 1].min(), rz_array[:, :, 1].max(), 64))
xyi = np.concatenate((xxi[:, :, np.newaxis], yyi[:, :, np.newaxis]), axis=2)
print np.shape(rz_array.reshape(64 * 64, 2))
print np.shape(xyi.reshape(64 * 64, 2))
print np.shape(frames[666, :, :].reshape(64 * 64))

# grid the data.
zi = griddata(rz_array.reshape(64 * 64, 2),
              frames[666, :, :].reshape(64 * 64),
              xyi.reshape(64 * 64, 2), method='linear')
# contour the gridded data, plotting dots at the randomly spaced data points.

print np.sum(np.isnan(zi))


plt.figure()
plt.contour(xxi[0, :], yyi[:, 0], zi.reshape(64, 64), 15, linewidths=0.5, colors='k')
plt.contourf(xxi[0, :], yyi[:, 0], zi.reshape(64, 64), 15, cmap=plt.cm.hot)
plt.colorbar()
plt.xlabel('R / cm')
plt.ylabel('z / cm')

plt.show()
