#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import skew, kurtosis

"""
Plot time series of individual pixels of the GPI diagnostic
"""

shotnr = 1120711005

frames_file = np.load('%d/%d_frames.npz' % ( shotnr, shotnr) )
frames  = frames_file['frames_normalized_mean']

frame0 = 15000

ix_list = np.array([50, 52, 54, 56, 58, 60])
jy_list = np.array([32, 32, 32, 32, 32, 32])
N = np.shape(frames)[0]

nx = np.size(ix_list)

plt.figure(0)
plt.title('# %d' % shotnr)
for ix, jy, idx in zip( ix_list, jy_list, np.arange( np.size(ix_list) ) ):
    timeseries = frames[frame0:, jy, ix]

    plt.figure(0)
    plt.plot( timeseries + float(idx)*3.0, label = 'z = %d' % jy )

    # Compute skewness and kurtosis
    pdf_s = skew( timeseries )
    pdf_k = kurtosis( timeseries )

    plt.figure(idx + 1)
    plt.title('(Z,R) = (%d,%d), S = %4f, F = %4f' % (ix, jy, pdf_s, pdf_k) )
    plt.hist( timeseries, bins = 100, log = True )

plt.show()