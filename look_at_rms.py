#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from phantom_helper import make_rz_array
from scipy.interpolate import griddata

"""
Compute and plot
1. <I>_RMS (R,Z)    RMS as a function of R,Z
2. <I>_t (R,Z)      Time avg. as a function of R,Z
3. S(I) (R,Z)       Skewness as a function of R,Z
4. K(I) (R,Z)       Kurtosis as a function of R,Z
"""

shotnr = 1120217008
datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr) )
frame_info = datafile['frame_info']
frame_range = datafile['framelist']
rz_array, transform_data = make_rz_array(frame_info)
#frame_range = np.arange(frame0, frame0+nframes)

rms  = datafile['frames_rms']
mean = datafile['frames_mean']


xxi, yyi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
xyi = np.concatenate( (xxi[:,:,np.newaxis], yyi[:,:,np.newaxis]), axis=2 )


print np.shape(rz_array)

# grid the data.
rms_grid  = griddata(rz_array.reshape(64*64, 2), rms.reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' )
mean_grid = griddata(rz_array.reshape(64*64, 2), mean.reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' )

plt.figure()
plt.title('#%d, $\\langle I \\rangle_t (R,z)$' % (shotnr) )
plt.xlabel('R / cm')
plt.ylabel('Z / cm')
plt.contour (xxi[0,:], yyi[:,0], mean_grid.reshape(64, 64), 15, linewidths=0.5, colors='k')
plt.contourf(xxi[0,:], yyi[:,0], mean_grid.reshape(64, 64), 15, cmap = plt.cm.hot)
plt.colorbar()


plt.figure()
plt.title('#%d, $\\langle I \\rangle_{RMS} (R,z)$' % ( shotnr ))
plt.xlabel('R / cm')
plt.ylabel('Z / cm')
plt.contour (xxi[0,:], yyi[:,0], rms_grid.reshape(64, 64), 15, linewidths=0.5, colors='k')
plt.contourf(xxi[0,:], yyi[:,0], rms_grid.reshape(64, 64), 15, cmap = plt.cm.hot)
plt.colorbar()


plt.figure()
plt.title('#%d, $\\langle I \\rangle_{t} / \\langle I \\rangle_{RMS} (R,z)$' % ( shotnr ))
plt.xlabel('R / cm')
plt.ylabel('Z / cm')
plt.contour (xxi[0,:], yyi[:,0], (mean_grid/rms_grid).reshape(64, 64), 15, linewidths=0.5, colors='k')
plt.contourf(xxi[0,:], yyi[:,0], (mean_grid/rms_grid).reshape(64, 64), 15, cmap = plt.cm.hot)
plt.colorbar()



plt.show()

#np.savez('%d/%d_rms.npz' % ( shotnr, shotnr ), rms=rms, shotnr = shotnr, transposed=True)


