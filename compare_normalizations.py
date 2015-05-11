#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
#import scipy.stats as s
import matplotlib.pyplot as plt

"""
Compare individual normalizations:

1.) I / <I>_t
2.) ( I - <I>_t ) / <I>_t
3.) ( I - <I>_t ) / <I>_RMS

"""

frames_file = np.load('1100803015/frames.npz')
frames_int = frames_file['frames_data']
frames_int[:,:,:] = frames_int[:,:,::-1]
# Convert to float datatype
frames = frames_int.astype('float')

t_idx = 23000

avg_t = np.sum( frames, axis=0 ) / float( np.shape(frames)[0] )
avg_rms = np.sum( (frames-avg_t)**2, axis=0 ) / np.sqrt( float( np.shape(frames)[0] ) )

plt.figure()
plt.subplot(221)
plt.title('Shot 1100803015, frame %d' % t_idx)
plt.contourf(frames[t_idx, :, :], cmap=plt.cm.hot, axes='equal')
plt.colorbar()

plt.subplot(222)
plt.title('$I / \\langle I \\rangle_t$')
plt.contourf(frames[t_idx, :, :] / avg_t, cmap = plt.cm.hot, axes='equal' )
plt.colorbar()

plt.subplot(223)
plt.title('$I - \\langle I \\rangle_t / \\langle I \\rangle_t$')
plt.contourf( ( frames[t_idx, :, :] - avg_t ) / avg_t , cmap = plt.cm.hot, axes='equal' )
plt.colorbar()

plt.subplot(224)
plt.title('$I - \\langle I \\rangle_t / \\langle I \\rangle_{RMS}$')
plt.contourf( ( frames[t_idx, :, :] - avg_t ) / avg_rms, cmap = plt.cm.hot, axes='equal'  )
plt.colorbar()

plt.show()



