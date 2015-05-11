#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cond_avg_2d import cond_avg_top_peak_surface


"""
Test conditional averaging routine
"""


frame0 = 30000
nframes = 50000
shotnr = 1100803006
#datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr) )
try:
    datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr), mmap_mode = 'c')
    frames = datafile['frames_normalized_mean']
    print 'Loaded frames for shot %d' % shotnr
except IOError:
    print 'Could not open file %d/%d_frames.npz' % (shotnr, shotnr)


px      = np.array([32,50])         # Choose a pixel well after the separatrix
minmax  = np.array([3., 6.])      # Peaks within 1.5 and 2.5 times above rms
lag     = 20
t       = np.arange(-lag, lag)
tseries = frames[:, px[0], px[1]]
time    = np.arange(0, nframes)

plt.figure()
plt.plot(time[:frame0], tseries[:frame0], 'r--')
plt.plot(time[frame0:], tseries[frame0:], 'k')

# Detect peaks
cavg_window, idx_events = cond_avg_top_peak_surface(frames, px, minmax, frame0, lag, rel_idx=True)
num_events = np.size(idx_events)

plt.figure()
for row in np.arange(0, np.shape(cavg_window)[0]):
    plt.plot(t, cavg_window[row,:])
    
plt.figure()
plt.title('#%d, Average blob form at %d,%d' % (shotnr, px[0], px[1]))
plt.errorbar(t, cavg_window.mean(axis=0), yerr=cavg_window.std(axis=0) )

plt.figure()
plt.subplot(121)
plt.title('Waiting time')
plt.xlabel('Frames between maxima')
plt.ylabel('Count')
plt.hist(np.abs(idx_events[1:] - idx_events[:-1]), bins=20)

plt.subplot(122)
plt.title('Blob amplitude')
plt.hist(tseries[idx_events], bins=20)
plt.xlabel('$I-\\langle I \\rangle_t / \\langle I \\rangle_{RMS}$ / a.u.')
plt.ylabel('Count')
print idx_events


# Plot blob realizations
for idx, event in enumerate(idx_events):
    fig = plt.figure(figsize=(10,5))
    fig.add_subplot(121)
    plt.title('#%d, pixel (%d,%d)' % (shotnr, px[0], px[1]))
    plt.xlabel('time')
    plt.ylabel('$I -\\langle I \\rangle_t / \\langle I \\rangle_{RMS} $')
    plt.plot(t + event, cavg_window[idx,:])

    ax = fig.add_subplot(122)
    cavg_surf = frames[ frame0 + event - lag : frame0 + event + lag,:,:].mean(axis=0)
    plt.title('average shape' )
    c = plt.contourf(cavg_surf, cmap=plt.cm.hot)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(c, cax=cax)

    plt.xlabel('x / px')
    plt.ylabel('y / px')
    
    F = plt.gcf()
    F.savefig('%d/%d_blob.eps' % (shotnr, event) )
    plt.close()
    
#    for i in np.arange(event - lag/5, events + lag/5):
#        plt.figure()
#        plt.title('#%d, t_idx=%d' % ( shotnr, i) )#dx_events[row] ) )
#        c = plt.contourf(frames[frame0+i, :, :], cmap=plt.cm.hot, axes='equal')
#        plt.plot(px[1], px[0], 'ko')
#        plt.colorbar()    
#        plt.xlabel('x / px')
#        plt.ylabel('y / px')
    

plt.show()

