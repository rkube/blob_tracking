#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

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
frame0 = 0000                          # Begin analysis at this frame
nframes = 1000                         # Number of frames to analyze
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

datafile = np.load('../test/test_frames_200.npz')
frames = datafile['frames']
frame_info = datafile['frame_info']

# Detect peaks
idx_events = detect_peak_3d(frames[frame0:frame0+nframes,:,:], trigger, minmax, 0, lag, rel_idx=False)
num_events = np.shape(idx_events)[0]
print '%d blob events detected' % ( num_events )

# Get R,z projection, grid data
rz_array, transform_data = make_rz_array(frame_info)
RRi, zzi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
Rzi = np.concatenate( (RRi[:,:,np.newaxis], zzi[:,:,np.newaxis]), axis=2 )
zi = griddata(rz_array.reshape(64*64, 2), frames[666,:,:].reshape( 64*64 ), Rzi.reshape( 64*64, 2 ), method='linear' )


for event in idx_events:
    try:
        I0 = event[0]
        t0 = event[1]
        R0 = event[2]
        z0 = event[3]
        event_frame = frames[t0,:,:]
        r0_last = R0
        z0_last = z0
        print 'Peak at (%d,%d): %f' % (R0, z0, frames[t0,R0,z0])
        xycom = np.zeros([2, tau_b + tau_a])
        rzcom = np.zeros([2, tau_b + tau_a])            # Rz com position at time t0 +- tau
        fwhm_Rz = np.zeros([2, tau_b + tau_a])
        amp = np.zeros([tau_b + tau_a])
        for t_idx, t in enumerate(np.arange(t0 - tau_b, t0 + tau_a)):
            event_frame = frames[t, :, :]
            labels = pm.label(event_frame > 0.6 * I0)
            blob_area = pm.blob( labels, 'area', output='data')
            blob_cent = pm.blob( labels, 'centroid', output='data')

            # In this frame, the region with centroid closes to the original blob is the new blob
            # Ignore regions that have less than 10% the orignial area
            min_dist = np.sqrt(64.*64.+64.*64.)
            for i in np.where( blob_area > 0.1 * max(blob_area) )[0]:
                # TODO: Vectorize loop
                dist = np.sqrt( (blob_cent[i,1]-z0_last)**2 + (blob_cent[i,0]-r0_last)**2 )
                if ( dist < min_dist ):
                    min_dist = dist
                    min_idx = i

            # Compute the x and y COM coordinates of the blob, store
            blob_mask = labels!=(min_idx+1)
            event_masked = np.ma.MaskedArray(labels, mask = blob_mask, fill_value=0)
            xycom[:,t_idx] = com(event_masked)
            rzcom[:,t_idx] = com_rz(event_masked, rz_array[:,:,0], rz_array[:,:,1])
#            print 'Region center at (%d,%d), com coordinates (x,y) = (%f,%f), (R,z) = (%f,%f)' %\
#                (blob_cent[i,0], blob_cent[i,1], xycom[0,t_idx], xycom[1,t_idx], rzcom[0,t_idx], rzcom[1,t_idx])

            # Compute FWHM of radial / poloidal cross section at COM coordinates
            xcom_off, ycom_off = np.round(xycom[:,t_idx]).astype('int')
            fwhm_rad_idx = fwhm(event_frame[ ycom_off, xcom_off - 16 : xcom_off + 16 ])
            fwhm_pol_idx = fwhm(event_frame[ ycom_off - 16 : ycom_off + 16, xcom_off ])
            fwhm_Rz[0,t_idx] = RRi[ ycom_off, xcom_off - 16 + fwhm_rad_idx[1][1] ] - RRi[ ycom_off, xcom_off - 16 + fwhm_rad_idx[1][0] ]
            fwhm_Rz[1,t_idx] = zzi[ ycom_off - 16 + fwhm_pol_idx[1][1] , xcom_off] - zzi[ ycom_off - 16 + fwhm_pol_idx[1][0] , xcom_off ]

#            print 'FWHM: (R,z) = (%f,%f)' % (fwhm_Rz[0,t_idx], fwhm_Rz[1,t_idx])
            amp[t_idx] = event_frame[r0_last, z0_last]
            r0_last, z0_last = blob_cent[min_idx, :]

#            fig = plt.figure()
#            fig.add_subplot(121, aspect='equal')
#            plt.title('frame %d' % t)
#            plt.contour (RRi[0,:], zzi[:,0], event_frame, 16, colors = 'k', linewidth=0.5 )
#            plt.contourf(RRi[0,:], zzi[:,0], event_frame, 16, cmap = plt.cm.hot )
#            plt.colorbar()
#            plt.plot(rzcom[0,:t_idx+1], rzcom[1,:t_idx+1], 'ko')
#            plt.xlabel('R / cm')
#            plt.ylabel('z / cm')
#
#            fig.add_subplot(122)
#            plt.title('Cross sections')
#            plt.plot(event_frame[ycom_off, xcom_off - 16 : xcom_off + 16 ], 'b-o', label='Radial cross, FWHM=%3.1f' % fwhm_Rz[0,t_idx])
#            plt.plot( fwhm_rad_idx[1], event_frame[ ycom_off, (fwhm_rad_idx[1] + xcom_off - 16).astype('int') ], 'b--' )
#            plt.plot(event_frame[ ycom_off - 16 : ycom_off + 16, xcom_off], 'g-o', label='Poloidal cross, FWHM=%3.1f' % fwhm_Rz[0,t_idx] )
#            plt.plot( fwhm_pol_idx[1], event_frame[ (fwhm_pol_idx[1] + ycom_off - 16).astype('int'), xcom_off ], 'g--' )
#            plt.legend(loc='upper right')
#
#            print 'New peak at (%d,%d)=%f' % (r0_last, z0_last, event_frame[r0_last, z0_last])

        print 'Amplitude: %f Radial velocity: %f m/s' % ( amp.mean(), 0.01*(rzcom[0,-1] - rzcom[0,0]) / float(tau_a+tau_b)/dt )
        # Plot poloidal vs. radial position
        blob_t = np.arange(t0-tau_b, t0+tau_a)
        plt.figure()
        plt.subplot(121)
        plt.plot(rzcom[0,:], rzcom[1,:])
        plt.xlabel('R / cm')
        plt.ylabel('z / cm')


        # Plot blob size as a function of time
        plt.subplot(122)
        plt.plot(blob_t, fwhm_Rz[0,:], label='rad')
        plt.plot(blob_t, fwhm_Rz[1,:], label='pol')
        plt.plot(blob_t, amp, label='Intensity')
        plt.xlabel('time / frames')
        plt.ylabel('size / cm')
        plt.legend()
    except:
        print 'Oooops'


plt.show()
