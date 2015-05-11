#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from misc.phantom_helper import make_rz_array
from scipy.stats import skew, kurtosis

"""
Plot radial skewness/kurtosis profiles from GPI data
"""

shotnr = 1100803015
frame0 = 15000

frames_file = np.load('%d/%d_frames.npz' % ( shotnr, shotnr) )
frames  = frames_file['frames_normalized_mean']

len = np.shape(frames)[0]

plt.figure()
plt.plot(np.arange(frame0), frames[:frame0, 32, 60], 'k')
plt.plot(np.arange(frame0, len), frames[frame0:, 32, 60], 'r')


rz_array, transform_data = make_rz_array(frames_file['frame_info'])
xxi, yyi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
xyi = np.concatenate( (xxi[:,:,np.newaxis], yyi[:,:,np.newaxis]), axis=2 )

pdf_s = skew( frames[frame0:, 32, :], axis = 0)
pdf_k = kurtosis( frames[frame0:, 32, :], axis = 0)

plt.figure()
plt.plot(xyi[32, :, 0], pdf_s, label = 'Skewness' )
plt.plot(xyi[32, :, 0], pdf_k, label = 'Kurtosis')
plt.legend(loc = 'upper left')
plt.xlabel('R / cm')
plt.ylabel('Skewness, Kurtosis')
plt.title('# %d, Single pixel statistics' % shotnr)

plt.show()