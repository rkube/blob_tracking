#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

from detect_peak import detect_peak_3d


"""
Test conditional averaging routine
"""

np.set_printoptions(linewidth=999999)
shotnr = 1100803006
try:
    datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr), mmap_mode = 'c')
    frames = datafile['frames_normalized_mean']
    print 'Loaded frames for shot %d' % shotnr
except IOError:
    print 'Could not open file %d/%d_frames.npz' % (shotnr, shotnr)

#datafile = np.load('../test/test_frames_200.npz')
#frames = datafile['frames']

frame0 = 20000                          # Begin analysis at this frame
nframes = 30000                         # Number of frames to analyze
minmax  = np.array([2.5, 4.5])          # Peaks within 1.5 and 2.5 times above rms
lag     = 20                            # Deadtime after a peak in which no blob is detected
time    = np.arange(0, nframes)         # Just the time, in frame numbers of the time series
trigger = np.array([40, 55, 8, 56])     # Triggerbox r_low, r_up, z_low, z_up
blobbox = np.array([8,8])               # Box used to determine blob size around single blob events
toffset = frame0 + lag                  # Total frame offset used in this script.


plt.figure()
plt.plot(frames[frame0:frame0+nframes, trigger[2]+8, trigger[0]+5])

# Detect peaks
idx_events = detect_peak_3d(frames[frame0:frame0+nframes,:,:], trigger, minmax, 0, lag, rel_idx=False)

for event in idx_events:
    print event

num_events = np.shape(idx_events)[0]
print '%d events' % num_events

print 'Removing events from the edge of the bounding box...'
idx_events = idx_events[idx_events['ridx'] > trigger[0] + 2]
idx_events = idx_events[idx_events['ridx'] < trigger[1] - 2]


#print 'returned thse peaks: ', idx_events
num_events = np.shape(idx_events)[0]
print '%d events' % num_events


plt.figure()
plt.hist( idx_events['value'], bins=50)
plt.title('Blob amplitude occurence')
plt.ylabel('Occurence')
plt.xlabel('Amplitude')

avg_blob = np.zeros([16,16])
Rr = np.arange(0,16)
zr = np.arange(0,16)

all_plot = False
for event in idx_events:
    t0 = event['tidx'] + frame0
    R0 = event['ridx'] # + trigger[0]
    z0 = event['zidx'] # + trigger[2]

    if all_plot:
        plt.figure()
        plt.subplot(221)
        plt.title('tidx=%d, ridx=%d zidx=%d' % (event['tidx']+frame0, event['ridx'], event['zidx']) )
        plt.contourf( frames[t0, :, : ], 32)
        plt.plot( R0 , z0 , 'ko')
     
        # Plot trigger box
        plt.plot( (trigger[1], trigger[0]), (trigger[2], trigger[2]), 'w')
        plt.plot( (trigger[1], trigger[0]), (trigger[3], trigger[3]), 'w')
        plt.plot( (trigger[0], trigger[0]), (trigger[2], trigger[3]), 'w')
        plt.plot( (trigger[1], trigger[1]), (trigger[2], trigger[3]), 'w')
    
        plt.colorbar()
        
    
        plt.subplot(222)
        plt.contourf(Rr, zr, frames[t0, z0 - 8 : z0+8, R0 -8 : R0 + 8], 32 )
        plt.colorbar()
    
        plt.subplot(223)
        plt.plot( frames[t0, z0, :], label='cut r' )
        plt.plot( frames[t0, :, R0], label='cut z')    
        plt.legend()

    avg_blob[:,:] = avg_blob[:,:] + frames[t0, z0 - 8 : z0+8, R0 -8 : R0 + 8]

avg_blob = avg_blob / float(num_events)

plt.figure()
plt.subplot(121)
plt.title('Average blob shape')
plt.contourf(avg_blob, 32)
plt.colorbar()

plt.subplot(122)
plt.plot(Rr, avg_blob[8,:], label='z cut')
plt.plot(zr, avg_blob[:,8], label='R cut')
plt.show()
