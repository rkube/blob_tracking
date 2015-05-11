#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt

"""
Plot time series of individual pixels of the gPI
"""

shotnr = 1120711005

frames_file = np.load('%d/%d_frames.npz' % ( shotnr, shotnr) )
frames  = frames_file['frames_normalized_mean']

frame0 = 15000

ix_list = np.array([32, 48, 56])
jy_list = np.array([32, 32, 32])
N = np.shape(frames)[0]

nx = np.size(ix_list)

plt.figure(0)
plt.title('# %d' % shotnr)
for ix, jy, idx in zip( ix_list, jy_list, np.arange( np.size(ix_list) ) ):
    plt.figure(0)
    plt.plot( frames[frame0:, jy, ix] + float(idx)*3.0, label = 'z = %d' % jy )

    plt.figure(idx + 1)
    plt.title('(Z,R) = (%d,%d)' % (ix, jy) )
    plt.hist( frames[frame0:, jy, ix], nbins = 100, log = True )

plt.show()