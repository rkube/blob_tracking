#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
from os.path import join

"""
Loads frames grabbed from MDS tree in memory
"""

def load_mdsframes(shotnr, path=None, fname='frames.npz', varname='frames_normalized_mean'):

    """
    Load phantom frames stored in a npz file.
    See /home/rkube/misc_scripts/save_phantom_frames.py

    input:
        shotnr:  int, Get phantom data from this show
        path:    string, Directory where the datafiles are stored
        varname: string, Variable to extract from the datafile

    output:
        frames:     ndarray, axis0: time, axis1: poloidal, axis2: radial
        frame_info: dictionary, Geometry of frames, see save_phantom_frames.py

    """

    # Set the default path
    if path is None:
        path = '/Users/ralph/source/blob_tracking/%10d' % (shotnr)
    if fname is None:
        fname = '%s/%10d_frames.npz' % (path, shotnr)

    fname = '%s/%s' % (path, fname)

    print 'Loading frames from %s/%d_frames.npz' % (path, shotnr)
    df = np.load(fname, mmap_mode='c')
    frames = df[varname]
    frame_info = df['frame_info']
    df.close()

    return frames, frame_info


def load_feltorframes(datadir, fname='frames_normalized.npz', varname='frames_normalized_rms'):
    """
    Load frames constructed from electron density of feltorSlab simulations.
    Script used to create datafiles: ~/uni/feltor_runs/make_frame_blobtracking.py
    """
    
    print 'Loading frames from %s/%s' % (datadir, fname)
    fname = join(datadir, fname)

    df = np.load(fname, mmap_mode='c')
    frames = df[varname]


    # Create a dummy geom object that can be passed to make_rz_array, find_sol_px, etc.
    geom = {'x_sol': df['x_sol'],
            'xrg': df['xrg_ip'],
            'yrg': df['yrg_ip'],
            'Nx': df['frames_normalized_rms'].shape[2],
            'Ny': df['frames_normalized_rms'].shape[1]}

    df.close()

    return frames, geom

#End of file load_mdsframes.py
