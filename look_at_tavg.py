#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt

frame0 = 0
nframes = 50000

datafile = np.load('1100803015/frames.npz')
frames = datafile['frames_data'][frame0:frame0+nframes, :, ::-1]
frame_range = np.arange(frame0, frame0+nframes)

print 'Computing <I>_t'

tavg = np.sum(frames, axis=0) / float(np.shape(frames)[0])
plt.figure()
plt.title('$\\langle I \\rangle_t (R,Z)$')
plt.contourf(tavg)
plt.colorbar()


plt.figure()
plt.title('$\\langle I \\rangle_t (R, Z = 32)$')
plt.plot(tavg[32,:])

plt.show()

