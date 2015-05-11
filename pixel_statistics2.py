#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
from plotting.figure_defs import set_rcparams_pres_small
import matplotlib as mpl
#mpl.use('AGG')
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt
from misc.phantom_helper import make_rz_array
from misc.load_mdsframes import load_mdsframes
from scipy.interpolate import griddata
from scipy.stats import kurtosis, skew
import idlsave

"""
Compute statistics for single pixels of the GPI array
"""


#shots = [ 1100803005, 1100803006, 1100803007, 1100803008, 1100803009, 1100803011, 1100803012, 1100803012,\
#    1100803015, 1100803020]
#shots = [ 1120217008, 1120217010, 1120217011, 1120217012, 1120217014, 1120217015, 1120217016, \
#    1120217017, 1120217018, 1120217019, 1120217020, 1120217021]
#shots = [ 1100811010 ]
shots = [ 1120711005, 1120711006, 1120711007, 1120711008, 1120711009, 1120711010, 1120711011 ]
#frame0 = 0    

for shotnr in shots:
    
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
    
    print np.shape(rz_array)
    
    r_flat = rz_array[:,:,0].reshape(64*64)
    z_flat = rz_array[:,:,1].reshape(64*64)
    
    # define grid.
    xxi, yyi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
    xyi = np.concatenate( (xxi[:,:,np.newaxis], yyi[:,:,np.newaxis]), axis=2 )
    # grid the data.
#    zi = griddata(rz_array.reshape(64*64, 2), frames[666,:,:].reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' )
    xxi_flat, yyi_flat = xxi.reshape(64*64), yyi.reshape(64*64)
    
    print 'shotnr %d' % shotnr
    print 'Separatrix at %f cm' % (s.rmid_sepx * 100.)
    print 'Limiter shadow at %f cm' % (s.rmid_lim * 100.)
    
    # Flat indices of points between separatrix and limiter shadow
    gap_idx = np.all( ((xxi_flat > s.rmid_sepx*100), (xxi_flat < s.rmid_lim*100)  ), axis=0)
    # 2d array indices of pixels between separatrix and limiter shadow
    gap_idx_unr = np.unravel_index(np.where(gap_idx)[0], (64,64) )
    gap_idx_unr = np.concatenate( (gap_idx_unr[0][:,np.newaxis], gap_idx_unr[1][:,np.newaxis]), axis=1)
    
    # Analyze only pixel at poloidal index 32
    n = np.size(np.where(gap_idx_unr[:,0] == 32)[0])
    # Take only 5 positions
    r_idx = np.where(gap_idx_unr[:,0] == 32)[0][ :: np.round(n / 5)+1 ]
    n_r = np.size( np.array(r_idx) )
    
#    # For frames, axis1 is the poloidal direction, axis2 is the radial direction
#    pixels = [frames[:, gap_idx_unr[r, 0], gap_idx_unr[r, 1] ] for r in r_idx]
#    # For rz_array, axis0 is the radial pixel, axis1 is the poloidal pixel
    radii  = [ rz_array[ gap_idx_unr[r, 0], gap_idx_unr[r, 1], 0 ] for r in r_idx]
    pol    = [ rz_array[ gap_idx_unr[r, 0], gap_idx_unr[r, 1], 1 ] for r in r_idx]
    
    
    # Take fluctuation and RMS from pre-computed array
    
    rms = np.array([ frames_rms[ gap_idx_unr[r, 0], gap_idx_unr[r, 1] ] for r in r_idx])
    mean  = np.array([ frames_mean [ gap_idx_unr[r, 0], gap_idx_unr[r, 1] ] for r in r_idx])    
    
    # Interpolate mean and rms on the new grid
    frames_mean_ip = griddata(rz_array.reshape(64*64, 2), frames_mean.reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' ).reshape(64,64)
    frames_rms_ip  = griddata(rz_array.reshape(64*64, 2), frames_rms.reshape( 64*64 ), xyi.reshape( 64*64, 2 ), method='linear' ).reshape(64,64)
    

    f = plt.figure( figsize = (15,8))
    f.text(0.5, 0.95, '#%d' % shotnr, ha = 'center')
    ax1 = f.add_subplot(2, 3, 1, aspect = 'equal')
    plt.title('$\\bar{I}$')
#    plt.contour ( xyi[:, :, 0], xyi[:, :, 1], frames_mean_ip, 32, linewidth = 0.5, colors = 'k')
    plt.contourf( xyi[: ,: ,0], xyi[:, :, 1], frames_mean_ip, 32, linewidth=0.5, cmap = plt.cm.hot)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')
    plt.colorbar()


    ax2 = f.add_subplot(2, 3, 2, aspect = 'equal')
    plt.title('$I_{\\mathrm{RMS}}$')
#    plt.contour ( xyi[:, :, 0], xyi[:, :, 1], frames_rms_ip, 32, linewidth = 0.5, colors = 'k')
    plt.contourf( xyi[: ,: ,0], xyi[:, :, 1], frames_rms_ip, 32, linewidth=0.5, cmap = plt.cm.hot)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')
    plt.colorbar()
    
    ax3 = f.add_subplot(2, 3, 3, aspect = 'equal')
    plt.title('$\\bar{I} / I_{\\mathrm{RMS}}$')
#    plt.contour ( xyi[:, :, 0], xyi[:, :, 1], (frames_rms_ip / frames_mean_ip), 32, linewidth = 0.5, colors = 'k')
    plt.contourf( xyi[: ,: ,0], xyi[:, :,1], (frames_rms_ip / frames_mean_ip), 32, linewidth=0.5, cmap = plt.cm.hot)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')
    plt.colorbar()
    
    
    f.add_subplot(2, 3, 4)
    plt.plot( xyi[32, :, 0], frames_mean[32,:] )
    plt.xlabel('R / cm')
    plt.ylabel('$\\bar{I}$')
    
    f.add_subplot(2, 3, 5)
    plt.plot( xyi[32, :, 0], frames_rms[32,:] )
    plt.xlabel('R / cm')
    plt.ylabel('$I_{RMS}$')
    
    f.add_subplot(2, 3, 6)
    plt.plot( xyi[32, :, 0], frames_rms[32,:] / frames_mean[32,:] )
    plt.xlabel('R / cm')
    plt.ylabel('$I_{RMS} / \\bar{I}$')
    
    f.savefig('%d/statistics/%d_fov_stats.eps' % (shotnr, shotnr))
    
plt.show()
