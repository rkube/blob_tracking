#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt

datafile = np.load('1100803015/frames.npz')
trange = np.arange(10000,50000)
frames = datafile['frames_data'][trange, :, ::-1]
mean16 = np.zeros([64, 64], dtype='int16')

print frames.dtype
mean = np.sum(frames, axis=0, dtype='uint64') / np.shape(frames)[0]

mean16[:,:] = mean[:,:]

plt.figure()
plt.title('1100803015, mean intensity')
plt.contourf(mean, cmap = plt.cm.hot )
plt.colorbar()

plt.show()

print mean16.dtype
np.savez('1100803015/frames_mean.npz', frames_mean = mean16 )

