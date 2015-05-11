#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt

"""
Plot time series of individual pixels.
"""

shotnr = 1100819016

frames_file = np.load('%d/%d_frames.npz' % ( shotnr, shotnr) )
frames_int = frames_file['frames_data']
frames_int[:,:,:] = frames_int[:,:,::-1]
# Convert to float
frames = frames_int.astype('float')

ix_list = [58, 58]
jy_list = [32, 31]
N = np.shape(frames_int)[0]

timeseries = np.zeros([40000, np.size(ix_list)])


for ix, jy, idx in zip( ix_list, jy_list, np.arange( np.size(ix_list) ) ):
    timeseries[:, idx] = frames[10000:, jy, ix]

plt.figure()
plt.subplot(211)
plt.title('Shot # %d' % shotnr)
plt.plot( timeseries[:, 0], label='%d, %d' % (ix_list[0], jy_list[0]) )
plt.plot( timeseries[:, 1], label='%d, %d' % (ix_list[1], jy_list[1]) )
plt.legend()

plt.subplot(212)
plt.hist( timeseries[:, 0], bins=100, log=True, histtype='step')
plt.hist( timeseries[:, 1], bins=100, log=True, histtype='step')


hist_i1, bins_i1 = np.histogram( timeseries[:,0], bins=100, normed=True )
hist_i2, bins_i2 = np.histogram( timeseries[:,1], bins=100, normed=True )


# Compute mean, rms
mean_1 = np.mean(timeseries[:,0])
mean_2 = np.mean(timeseries[:,1])

rms_1 = np.sum( timeseries[:,0] * timeseries[:,0] / float(N) )
rms_2 = np.sum( timeseries[:,1] * timeseries[:,1] / float(N) )

# Plot normalized time series 
plt.figure()
plt.subplot(211)
plt.title('$ I / \\langle I \\rangle_t$')
plt.plot(timeseries[:,0] / mean_1)
plt.plot(timeseries[:,1] / mean_2)
plt.subplot(212)
# Compute normalized histograms
print np.shape(bins_i1), np.shape(hist_i1)
plt.semilogy( bins_i1[1:] / mean_1, hist_i1 * mean_1 )
plt.semilogy( bins_i2[1:] / mean_2, hist_i2 * mean_2 )
#h1, bins1 = np.histogram( timeseries[:,0] / mean_1, bins=100, density = True )
#h2, bins2 = np.histogram( timeseries[:,1] / mean_2, bins=100, density = True )
#plt.semilogy(bins1[1:] , hist_i1 * mean_1 )
#plt.semilogy(bins2[1:] , hist_i2 * mean_2 )
sum1 = np.sum( hist_i1 * (bins_i1[1:] - bins_i1[:-1]) )
sum2 = np.sum( hist_i2 * (bins_i2[1:] - bins_i2[:-1]) )
print 'Normalization 1: Int(PDF(1)) = %f, Int(PDF(2)) = %f' % (sum1, sum2)

plt.figure()
plt.subplot(211)
plt.title('$ I - \\langle I \\rangle_t / \\langle I \\rangle_{t}$')
plt.plot( (timeseries[:,0] - mean_1 ) / mean_1 , label='Mean = %f' % mean_1 )
plt.plot( (timeseries[:,1] - mean_2 ) / mean_2 , label='Mean = %f' % mean_2 )
plt.legend()
plt.subplot(212)
plt.semilogy( (bins_i1[1:] - mean_1) / mean_1, hist_i1 * mean_1 )
plt.semilogy( (bins_i2[1:] - mean_2) / mean_2, hist_i2 * mean_2 )
#h1, bins1 = np.histogram( (timeseries[:,0] - mean_1 ) / mean_1 , bins = 100, density = True ) 
#h2, bins2 = np.histogram( (timeseries[:,1] - mean_2 ) / mean_2 , bins = 100, density = True) 
#plt.semilogy( bins1[1:], h1 )
#plt.semilogy( bins2[1:], h2 )
sum1 = np.sum( hist_i1 * (bins_i1[1:] - bins_i1[:-1]) )
sum2 = np.sum( hist_i2 * (bins_i2[1:] - bins_i2[:-1]) )
print 'Normalization 2: Int(PDF(1)) = %f, Int(PDF(2)) = %f' % (sum1, sum2)


plt.figure()
plt.subplot(211)
plt.title('$ I - \\langle I \\rangle_t / \\langle I \\rangle_{RMS}$')
plt.plot( (timeseries[:,0] - mean_1 ) / rms_1 , label='RMS = %f' % rms_1 )
plt.plot( (timeseries[:,1] - mean_2 ) / rms_2 , label='RMS = %f' % rms_2 )
plt.legend()
plt.subplot(212)
plt.semilogy( ( bins_i1[1:] - mean_1 ) / rms_1, hist_i1 * rms_1 )
plt.semilogy( ( bins_i2[1:] - mean_2 ) / rms_2, hist_i2 * rms_2 )
#h1, bins1 = np.histogram(  (timeseries[:,0] - mean_1 ) / rms_1, bins = 100, density = True )
#h2, bins2 = np.histogram(  (timeseries[:,1] - mean_2 ) / rms_2, bins = 100, density = True )
#plt.semilogy( bins1[1:], h1 )
#plt.semilogy( bins2[1:], h2 )
sum1 = np.sum( hist_i1 * (bins_i1[1:] - bins_i1[:-1]) )
sum2 = np.sum( hist_i2 * (bins_i2[1:] - bins_i2[:-1]) )
print 'Normalization 3: Int(PDF(1)) = %f, Int(PDF(2)) = %f' % (sum1, sum2)


plt.show()

