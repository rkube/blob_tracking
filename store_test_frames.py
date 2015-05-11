#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from sys import argv

"""
Create a test data set for a given shot
"""

shotnr = 1111208021
nframes = 1000
frame0 = 8000
scriptname = argv[0]

datadir = '/Volumes/My Book Thunderbolt Duo/cmod_data/phantom/data/'
outputdir = '/Volumes/My Book Thunderbolt Duo/cmod_data/phantom/data/'

df_infname = '%s/%10d_frames_normalized.npz' % (datadir, shotnr)
df_outfname = '%s/%10d_testframes.npz' % (datadir, shotnr)

try:
    #datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr), mmap_mode = 'c')
    df = np.load(df_infname, mmap_mode='r')
    frames = df['frames_normalized_rms']
    frame_info = df['frame_info']
    print 'Loaded frames for shot %d' % shotnr
    df.close()
except IOError:
    print 'Could not open file %s' % (df_infname)

R = 50
Z = 32

test_frames = frames[frame0:frame0 + nframes, :, :]
t = np.arange(0, np.shape(frames)[0])
plt.figure()

plt.title('#%10d at pixel (R,Z) = (%d,%d)' % (shotnr, R, Z))
plt.plot(t[:frame0], frames[:frame0, Z, R], 'k')
plt.plot(t[frame0:frame0 + nframes],
         frames[frame0:frame0 + nframes, Z, R], 'r')
plt.plot(t[frame0 + nframes:], frames[frame0 + nframes:, Z, R], 'k')

np.savez(df_outfname, frames=test_frames, shotnr=shotnr,
         nframes=nframes, frame0=frame0, frame_info=frame_info,
         scriptname=scriptname)

print '#%10d: stored test frames in %s' % (shotnr, df_outfname)

plt.show()

#End of file store_test_frames.py
