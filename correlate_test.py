#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from misc.correlate import correlate
from sys import exit


"""
Correlate timeseries of two pixels. Trigger on large amplitudes.
Try estimating radial blob velocity from this.

"""

shotnr = 1120711010
frame0 = 20000
nframes = 25000
wlen = 10           # Window length for correlation analysis
gpi_fps = 390.8e3  # Phantom camera frame rate
tau = 1. / gpi_fps
dx = 0.061 / 64.
tau_max = 10

num_blobs = 100
blob_vel = np.zeros(num_blobs)

frames_file = np.load('%d/%d_frames.npz' % ( shotnr, shotnr) )
frames = frames_file['frames_normalized_mean']
ts_pixel1 = frames[ :, 48, 53 ]
ts_pixel2 = frames[ :, 48, 55 ]
ts_pixel3 = frames[ :, 48, 57 ]
ts_pixel4 = frames[ :, 48, 59 ]
ts_pixel5 = frames[ :, 48, 61 ]

#np.savez('corr_ts.npz', ts_pixel1 = ts_pixel1, ts_pixel2 = ts_pixel2, ts_pixel3 = ts_pixel3, ts_pixel4 = ts_pixel4, ts_pixel5 = ts_pixel5 )
#df = np.load('test/corr_ts.npz')
#ts_pixel1 = df['ts_pixel1']
#ts_pixel2 = df['ts_pixel2']
#ts_pixel3 = df['ts_pixel3']
#ts_pixel4 = df['ts_pixel4']
#ts_pixel5 = df['ts_pixel5']

plt.figure()
plt.plot( np.arange( frame0 ), ts_pixel1[:frame0], 'k' )
plt.plot( np.arange( frame0, frame0 + nframes ), ts_pixel1[frame0:frame0 + nframes], 'r' )
plt.plot( np.arange( frame0 + nframes, np.size(ts_pixel1)), ts_pixel1[frame0 + nframes:], 'k' )

plt.plot( np.arange( frame0 ), ts_pixel2[:frame0] + 3.0, 'k' )
plt.plot( np.arange( frame0, frame0 + nframes ), ts_pixel2[frame0:frame0 + nframes] + 3.0, 'r')
plt.plot( np.arange( frame0 + nframes, np.size(ts_pixel1)), ts_pixel2[frame0 + nframes:] + 3.0, 'k' )

plt.plot( np.arange( frame0 ), ts_pixel3[:frame0] + 6.0, 'k' )
plt.plot( np.arange( frame0, frame0 + nframes ), ts_pixel3[frame0:frame0 + nframes] + 6.0, 'r' )
plt.plot( np.arange( frame0 + nframes, np.size(ts_pixel1)), ts_pixel3[frame0 + nframes:] + 6.0, 'k' )

plt.plot( np.arange( frame0 ), ts_pixel4[:frame0] + 9.0, 'k' )
plt.plot( np.arange( frame0, frame0 + nframes ), ts_pixel4[frame0:frame0 + nframes]+ 6.0, 'r' )
plt.plot( np.arange( frame0 + nframes, np.size(ts_pixel1)), ts_pixel4[frame0 + nframes:] + 6.0, 'k' )

plt.plot( np.arange( frame0 ), ts_pixel5[:frame0] + 9.0, 'k' )
plt.plot( np.arange( frame0, frame0 + nframes ), ts_pixel5[frame0:frame0 + nframes] + 9.0, 'r' )
plt.plot( np.arange( frame0 + nframes, np.size(ts_pixel1)), ts_pixel5[frame0 + nframes:] + 9.0, 'k' )

#plt.show()


ts_pixel1 = ts_pixel1[frame0 : frame0 + nframes]
ts_pixel2 = ts_pixel2[frame0 : frame0 + nframes]
ts_pixel3 = ts_pixel3[frame0 : frame0 + nframes]
ts_pixel4 = ts_pixel4[frame0 : frame0 + nframes]
ts_pixel5 = ts_pixel5[frame0 : frame0 + nframes]

# Take the 100 largest blobs and estimate their velocity
ts1_sortidx = ts_pixel1.argsort()[-num_blobs:]
plt.figure()
plt.plot(ts_pixel1[ts1_sortidx])


for idx, max_idx in enumerate(ts1_sortidx):
    if ( max_idx == -1 ):
        print 'Index was blanked out previously, skipping to next index'
        continue
    elif ( max_idx < wlen ):
        print 'Too close too boundaries for full correlation, skipping to next index'
        continue
    # Blank out all other peaks occuring within +- 10 frames
    print 'before:', max_idx, ts1_sortidx
    close_peak_indices = np.squeeze(np.argwhere( np.abs(ts1_sortidx - max_idx) < 10 ))
    print 'close_peak_indices:', close_peak_indices, ' entries:', ts1_sortidx[ close_peak_indices ] 
    ts1_sortidx[ close_peak_indices ] = -1
    print 'after:', max_idx, ts1_sortidx    

    print max_idx
    fig = plt.figure()
    plt.subplot(211)
    plt.title('max_idx = %d' % max_idx)
    plt.xlabel('frame no.')
    plt.ylabel('I tilde')
    plt.plot( np.arange( frame0 + max_idx - wlen, frame0 + max_idx + wlen), ts_pixel1[ max_idx - wlen : max_idx + wlen ] )
    plt.plot( np.arange( frame0 + max_idx - wlen, frame0 + max_idx + wlen), ts_pixel2[ max_idx - wlen : max_idx + wlen ] )
    plt.plot( np.arange( frame0 + max_idx - wlen, frame0 + max_idx + wlen), ts_pixel3[ max_idx - wlen : max_idx + wlen ] )
    plt.plot( np.arange( frame0 + max_idx - wlen, frame0 + max_idx + wlen), ts_pixel4[ max_idx - wlen : max_idx + wlen ] )
    plt.plot( np.arange( frame0 + max_idx - wlen, frame0 + max_idx + wlen), ts_pixel5[ max_idx - wlen : max_idx + wlen ] )
   
    
    plt.subplot(212)
    plt.xlabel('Time lag tau')
    plt.ylabel('Correlation amplitude')

    tau_range = np.arange( -tau_max, tau_max )

    # Compute the correlation between the timeseries of neighbouring pixels. The maximum of
    # the correlation amplitude is used to compute the radial blob velocity.
    # Limit the neighbourhood in which the peak correlation amplitude may be to +- 10 frames.
    c11 = correlate( ts_pixel1[ max_idx - wlen - 1: max_idx + wlen + 1], ts_pixel1[ max_idx - wlen - 1 : max_idx + wlen + 1], 2*wlen )
    c11 = c11[ 2*wlen - tau_max : 2*wlen + tau_max]
    plt.plot(tau_range, c11)
    plt.plot(tau_range[c11.argmax()], c11.max(), 'ko')
    
    c12 = correlate( ts_pixel1[ max_idx - wlen - 1: max_idx + wlen + 1], ts_pixel2[ max_idx - wlen - 1 : max_idx + wlen + 1], 2*wlen )
    c12 = c12[ 2*wlen - tau_max : 2*wlen + tau_max]
    max_c12 = c12.argmax()
    plt.plot(tau_range, c12)
    plt.plot(tau_range[c12.argmax()], c12.max(), 'ko')
        
    c13 = correlate( ts_pixel1[ max_idx - wlen - 1: max_idx + wlen + 1], ts_pixel3[ max_idx - wlen - 1: max_idx + wlen +1 ], 2*wlen )
    c13 = c13[ 2*wlen - tau_max : 2*wlen + tau_max]
    max_c13 = c13.argmax()
    plt.plot(tau_range, c13)
    plt.plot(tau_range[c13.argmax()], c13.max(), 'ko')
    
    c14 = correlate( ts_pixel1[ max_idx - wlen - 1: max_idx + wlen + 1], ts_pixel3[ max_idx - wlen - 1: max_idx + wlen + 1], 2*wlen )
    c14 = c14[ 2*wlen - tau_max : 2*wlen + tau_max]
    max_c14 = c14.argmax()
    plt.plot(tau_range, c14)
    plt.plot(tau_range[c14.argmax()], c14.max(), 'ko')
    
    c15 = correlate( ts_pixel1[ max_idx - wlen - 1: max_idx + wlen + 1], ts_pixel5[ max_idx - wlen - 1: max_idx + wlen + 1], 2*wlen )
    c15 = c15[ 2*wlen - tau_max : 2*wlen + tau_max]
    max_c15 = c15.argmax()
    plt.plot(tau_range, c15)
    plt.plot(tau_range[c15.argmax()], c15.max(), 'ko')

    fig.savefig('%d/vrad_correlation/%d_frame%05d.png' % ( shotnr, shotnr, max_idx ) )
    plt.close()

   
    # Estimate radial blob velocity by propagation of correlation amplitude
    v_c12 = gpi_fps * 2.0*dx / (2*wlen - max_c12)
    v_c13 = gpi_fps * 4.0*dx / (2*wlen - max_c13)
    v_c14 = gpi_fps * 6.0*dx / (2*wlen - max_c14)
    v_c15 = gpi_fps * 8.0*dx / (2*wlen - max_c15)
    
    print 'Blob velocities from correlation method:'
    print 'px1 - px2: %f, px1 - px3: %f, px1 - px4: %f, px1 - px5: %f' % (v_c12, v_c13, v_c14, v_c15 )
    print 'mean: %f' % np.mean( np.array([v_c12, v_c13, v_c14, v_c15]) )
    blob_vel[idx] = np.mean( np.array([v_c12, v_c13, v_c14, v_c15]) )

blob_vel = blob_vel[blob_vel != 0]
plt.figure()
plt.plot(blob_vel, '.')
plt.xlabel('Blob event no.')
plt.ylabel('Radial velocity m/s')

print '================================================================================='
print 'mean over all blobs: %f' % blob_vel.mean()

plt.show()
