#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
from figure_defs import set_rcparams_pres_small
import matplotlib as mpl
#mpl.use('AGG')
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt
from phantom_helper import make_rz_array
from load_mdsframes import load_mdsframes
from scipy.interpolate import griddata
import idlsave

"""
Compute statistics for single pixels of the GPI array
"""


#shots = [ 1100803005, 1100803006, 1100803007, 1100803008, 1100803009, 1100803011, 1100803012, 1100803012,\
#    1100803015, 1100803020]

shots = [ 1100803005, 1100803015 ]
label_str = ['$n_e / n_G = 0.15$', '$n_e / n_G = 0.31$']

f = plt.figure( figsize = (6,6))
for idx, shotnr in enumerate(shots):
    
    datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr), mmap_mode = 'c')
    frame_info  = datafile['frame_info']
    frames_mean = datafile['frames_mean']
    frames_rms  = datafile['frames_rms']

    try:
        s = idlsave.read('%d/%d_separatrix.sav' % (shotnr, shotnr) , verbose = False)
    except:
        print 'Loading default separatrix values for #1120217008'
        s = idlsave.read('1120217008/1120217008_separatrix.sav', verbose = False)
    
    rz_array, transform_data = make_rz_array(frame_info)

    # define grid.
    xxi, yyi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
    xyi = np.concatenate( (xxi[:,:,np.newaxis], yyi[:,:,np.newaxis]), axis=2 )
    
    print 'shotnr %d' % shotnr
    print 'Separatrix at %f cm' % (s.rmid_sepx * 100.)
    print 'Limiter shadow at %f cm' % (s.rmid_lim * 100.)
    
    # Interpolate mean and rms on the new grid
    frames_mean_ip = griddata(rz_array.reshape(64*64, 2), frames_mean.reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' ).reshape(64,64)
    frames_rms_ip  = griddata(rz_array.reshape(64*64, 2), frames_rms.reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' ).reshape(64,64)
    
    plt.plot( xyi[32, :, 0], frames_rms[32,:] / frames_mean[32,:], label= label_str[idx] )

    plt.xlabel('R / cm', fontsize='24')
    plt.ylabel('$I_{RMS} / \\bar{I}$', fontsize='24')
    
plt.axvspan(s.rmid_sepx*100, s.rmid_lim*100, facecolor='0.5', alpha=0.5)
plt.xlim((88, 92))
plt.ylim((0.1, 0.55))
leg = plt.legend()
leg.get_frame().set_alpha(0.4)
    
plt.show()
