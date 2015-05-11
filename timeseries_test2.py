#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt


"""
Similar timeseries_test.py

Reads timeseries from two CMOD shots from the same pixel and compares the PDF's

"""

# Read timeseries from GPI frame data and store them together with the
# defined parameters below in a npz file. Once the lines below are run 
# for a given parameter set, read the data from the npz file to save tim.
#
shotnrs = np.array([1100803008, 1100803009, 1100803007])
ix, jy = 50, 32     # Position of pixel where we take timeseries
N = 50000           # Number of elements in timeseries
start = 10000       # Cutoff for timeseries
timeseries = np.zeros([N-start, 3])

for idx, shotnr in enumerate(shotnrs):
    print 'Reading file %d/%d_frames.npz' % (shotnr, shotnr)
    frames_file = np.load('%d/%d_frames.npz' % (shotnr, shotnr) )
    print 'Transposed = ', frames_file['transposed']
    if ( frames_file['transposed'] == False ):
        timeseries[:, idx] = frames_file['frames_data'][start:,jy,64-ix].astype('float')
    else:
        timeseries[:, idx] = frames_file['frames_data'][start:,jy,ix].astype('float')

# End data input. Store to npz file
np.savez('script_data/timeseries_test2.npz', timeseries = timeseries, ix=ix, jy=jy, N=N, start=start, shotnrs=shotnrs) 

# Read only the timeseries defined above from the npz file
#
#datafile = np.load('script_data/timeseries_test2.npz')
#shotnrs = datafile['shotnrs']
#ix, jy = datafile['ix'], datafile['jy']
#N = datafile['N']
#start = datafile['start']
#timeseries = datafile['timeseries']


# Plot time series and their pdf
plt.figure()
plt.subplot(421)
plt.plot( timeseries[:,0], label='Shot %d' % shotnrs[0] )
plt.plot( timeseries[:,1], label='Shot %d' % shotnrs[1] )
plt.plot( timeseries[:,2], label='Shot %d' % shotnrs[2] )
plt.legend()

plt.subplot(422)
plt.hist( timeseries[:,0], bins = 100, log=True, histtype='step')
plt.hist( timeseries[:,1], bins = 100, log=True, histtype='step')
plt.hist( timeseries[:,2], bins = 100, log=True, histtype='step')

# Make historgrams
hist_i1, bins_i1 = np.histogram( timeseries[:,0], bins=100, normed=True )
hist_i2, bins_i2 = np.histogram( timeseries[:,1], bins=100, normed=True )
hist_i3, bins_i3 = np.histogram( timeseries[:,2], bins=100, normed=True )
# Compute mean and PDF
mean_1 = np.mean(timeseries[:,0])
mean_2 = np.mean(timeseries[:,1])
mean_3 = np.mean(timeseries[:,2])

rms_1 = np.sum( timeseries[:,0] * timeseries[:,0] / float(N-start) )
rms_2 = np.sum( timeseries[:,1] * timeseries[:,1] / float(N-start) )
rms_3 = np.sum( timeseries[:,2] * timeseries[:,2] / float(N-start) )


# Plot normalized time series 
#plt.figure()
plt.subplot(423)
plt.title('$ I / \\langle I \\rangle_t$')
plt.plot(timeseries[:,0] / mean_1)
plt.plot(timeseries[:,1] / mean_2)
plt.plot(timeseries[:,2] / mean_3)
plt.subplot(424)
plt.semilogy( bins_i1[1:] / mean_1, hist_i1 * mean_1 )
plt.semilogy( bins_i2[1:] / mean_2, hist_i2 * mean_2 )
plt.semilogy( bins_i3[1:] / mean_3, hist_i3 * mean_3 )
sum1 = np.sum( hist_i1 * (bins_i1[1:] - bins_i1[:-1]) )
sum2 = np.sum( hist_i2 * (bins_i2[1:] - bins_i2[:-1]) )
sum3 = np.sum( hist_i3 * (bins_i3[1:] - bins_i3[:-1]) )
print 'Normalization 1: Int(PDF(1)) = %f, Int(PDF(2)) = %f' % (sum1, sum2)

# Plot timeseries normalized with RMS
#plt.figure()
plt.subplot(425)
plt.title('$ I - \\langle I \\rangle_t / \\langle I \\rangle_{t}$')
plt.plot( (timeseries[:,0] - mean_1 ) / mean_1 , label='Mean = %f' % mean_1 )
plt.plot( (timeseries[:,1] - mean_2 ) / mean_2 , label='Mean = %f' % mean_2 )
plt.plot( (timeseries[:,1] - mean_3 ) / mean_3 , label='Mean = %f' % mean_3 )
plt.legend()
plt.subplot(426)
plt.semilogy( (bins_i1[1:] - mean_1) / mean_1, hist_i1 * mean_1 )
plt.semilogy( (bins_i2[1:] - mean_2) / mean_2, hist_i2 * mean_2 )
plt.semilogy( (bins_i3[1:] - mean_3) / mean_3, hist_i3 * mean_3 )
sum1 = np.sum( hist_i1 * (bins_i1[1:] - bins_i1[:-1]) )
sum2 = np.sum( hist_i2 * (bins_i2[1:] - bins_i2[:-1]) )
print 'Normalization 2: Int(PDF(1)) = %f, Int(PDF(2)) = %f' % (sum1, sum2)


#plt.figure()
plt.subplot(427)
plt.title('$ I - \\langle I \\rangle_t / \\langle I \\rangle_{RMS}$')
plt.plot( (timeseries[:,0] - mean_1 ) / rms_1 , label='RMS = %f' % rms_1 )
plt.plot( (timeseries[:,1] - mean_2 ) / rms_2 , label='RMS = %f' % rms_2 )
plt.plot( (timeseries[:,1] - mean_3 ) / rms_3 , label='RMS = %f' % rms_3 )
plt.legend()
plt.subplot(428)
plt.semilogy( ( bins_i1[1:] - mean_1 ) / rms_1, hist_i1 * rms_1 )
plt.semilogy( ( bins_i2[1:] - mean_2 ) / rms_2, hist_i2 * rms_2 )
plt.semilogy( ( bins_i3[1:] - mean_3 ) / rms_3, hist_i3 * rms_3 )
sum1 = np.sum( hist_i1 * (bins_i1[1:] - bins_i1[:-1]) )
sum2 = np.sum( hist_i2 * (bins_i2[1:] - bins_i2[:-1]) )
print 'Normalization 3: Int(PDF(1)) = %f, Int(PDF(2)) = %f' % (sum1, sum2)



plt.show()

