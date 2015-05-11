#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

from IPython import embed
import numpy as np
import pymorph as pm
import matplotlib.pyplot as plt
from detect_peak import detect_peak_3d
from helper_functions import com, com_rz, fwhm
from phantom_helper import make_rz_array
from scipy.interpolate import griddata


"""
Check out how image segmenting works
"""
np.set_printoptions(linewidth=999999)
shotnr  = 1100803015
frame0 = 20000                          # Begin analysis at this frame. 
nframes = 30000                         # Number of frames to analyze
minmax  = np.array([2.5, 10.0])         # Peaks within 1.5 and 2.5 times above rms
lag     = 20                            # Deadtime after a peak in which no blob is detected
time_arr= np.arange(0, nframes)         # Just the time, in frame numbers of the time series
trigger = np.array([40, 50, 10, 53])     # Triggerbox r_low, r_up, z_low, z_up
blobbox = np.array([8,8])               # Box used to determine blob size around single blob events
toffset = frame0 + lag                  # Total frame offset used in this script.
tau_b   = 2                             # Frames to plot before/after blob enters trigger box
tau_a   = 5
nbins   = 10                            # Bins for amplitude sorting
# 1 frame is 2.5mus
dt = 1./400000.

try:
    datafile = np.load('../blob_tracking/%d/%d_frames.npz' % (shotnr, shotnr), mmap_mode = 'c')
    frames = datafile['frames_normalized_mean']
    print 'Loaded frames for shot %d' % shotnr
except IOError:
    print 'Could not open file %d/%d_frames.npz' % (shotnr, shotnr)
    datafile = np.load('../test/test_frames_200.npz')
    frames = datafile['frames']
frame_info = datafile['frame_info']

embed()