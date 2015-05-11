#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-


import numpy as np
import pymorph as pm
import matplotlib.pyplot as plt
from detect_peak import detect_peak_3d
from helper_functions import com, com_rz, fwhm
from phantom_helper import make_rz_array
from scipy.interpolate import griddata

"""
Detect peaks
The detect_peak_3d returns time indices of blob events relative for the
array passed to it. Remember to add frame0 to idx_event[t0] to
translate to the frame indexing used in this script.
"""

frames = np.zeros([50000, 64, 64])
frame0 = 10000
nframes = 20000
idx_events = detect_peak_3d(frames[frame0:frame0 + nframes, :, :], trigger,
                            minmax, 0, lag, rel_idx=False)
num_events = np.shape(idx_events)[0]
event_ctr = np.ones([num_events])
print '%d blob events detected' % ( num_events )

# Get R,z projection, grid data
rz_array, transform_data = make_rz_array(frame_info)
RRi, zzi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
Rzi = np.concatenate( (RRi[:,:,np.newaxis], zzi[:,:,np.newaxis]), axis=2 )
zi = griddata(rz_array.reshape(64*64, 2), frames[666,:,:].reshape( 64*64 ), Rzi.reshape( 64*64, 2 ), method='linear' )


# Average amplitude, velocity, and length for each blob event
blob_amp = np.zeros([num_events])
blob_vel = np.zeros([num_events])
blob_lp  = np.zeros([num_events])
blob_lr  = np.zeros([num_events])

for idx, event in enumerate(idx_events[:2]):
#    try:
    for bob in np.arange(0,1):
        I0 = event[0]
        t0 = event[1] + frame0
        R0 = event[2]
        z0 = event[3]
        event_frame = frames[t0,:,:]
        r0_last = R0
        z0_last = z0
        print 'Frame %d, peak at (%d,%d): %f or %f' % (t0, R0, z0, frames[t0,R0,z0], I0)
        xycom = np.zeros([2, tau_b + tau_a])
        rzcom = np.zeros([2, tau_b + tau_a])            # Rz com position at time t0 +- tau
        fwhm_Rz = np.zeros([2, tau_b + tau_a])
        amp = np.zeros([tau_b + tau_a])

        plt.figure()
        plt.contourf( frames[t0, :, :], cmap = plt.cm.hot )
        plt.colorbar()
        plt.plot( z0, R0, 'ko' )

        for t_idx, t in enumerate(np.arange(t0 - tau_b, t0 + tau_a)):
            event_frame = frames[t, :, :]
            labels = pm.label(event_frame > 0.6 * I0)
            blob_area = pm.blob( labels, 'area', output='data')
            blob_cent = pm.blob( labels, 'centroid', output='data')

            print blob_area, blob_cent
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
            print 'Region center at (%d,%d), com coordinates (x,y) = (%f,%f), (R,z) = (%f,%f)' %\
                (blob_cent[i,0], blob_cent[i,1], xycom[0,t_idx], xycom[1,t_idx], rzcom[0,t_idx], rzcom[1,t_idx])

            # Compute FWHM of radial / poloidal cross section at COM coordinates
            xcom_off, ycom_off = np.round(xycom[:,t_idx]).astype('int')
            fwhm_rad_idx = fwhm(event_frame[ ycom_off, xcom_off - 16 : xcom_off + 16 ])
            fwhm_pol_idx = fwhm(event_frame[ ycom_off - 16 : ycom_off + 16, xcom_off ])
            fwhm_Rz[0,t_idx] = RRi[ ycom_off, xcom_off - 16 + fwhm_rad_idx[1][1] ] - RRi[ ycom_off, xcom_off - 16 + fwhm_rad_idx[1][0] ]
            fwhm_Rz[1,t_idx] = zzi[ ycom_off - 16 + fwhm_pol_idx[1][1] , xcom_off] - zzi[ ycom_off - 16 + fwhm_pol_idx[1][0] , xcom_off ]

            print 'FWHM: (R,z) = (%f,%f)' % (fwhm_Rz[0,t_idx], fwhm_Rz[1,t_idx])
            amp[t_idx] = event_frame[r0_last, z0_last]
            r0_last, z0_last = blob_cent[min_idx, :]

            fig = plt.figure()
            fig.add_subplot(121, aspect='equal')
            plt.title('frame %d' % t)
            plt.contour (RRi[0,:], zzi[:,0], event_frame, 16, colors = 'k', linewidth=0.5 )
            plt.contourf(RRi[0,:], zzi[:,0], event_frame, 16, cmap = plt.cm.hot )
            plt.colorbar()
            plt.plot(rzcom[0,:t_idx+1], rzcom[1,:t_idx+1], 'ko')
            plt.xlabel('R / cm')
            plt.ylabel('z / cm')

            fig.add_subplot(122)
            plt.title('Cross sections')
            plt.plot(event_frame[ycom_off, xcom_off - 16 : xcom_off + 16 ], 'b-o', label='Radial cross, FWHM=%3.1f' % fwhm_Rz[0,t_idx])
            plt.plot( fwhm_rad_idx[1], event_frame[ ycom_off, (fwhm_rad_idx[1] + xcom_off - 16).astype('int') ], 'b--' )
            plt.plot(event_frame[ ycom_off - 16 : ycom_off + 16, xcom_off], 'g-o', label='Poloidal cross, FWHM=%3.1f' % fwhm_Rz[0,t_idx] )
            plt.plot( fwhm_pol_idx[1], event_frame[ (fwhm_pol_idx[1] + ycom_off - 16).astype('int'), xcom_off ], 'g--' )
            plt.legend(loc='upper right')

            print 'New peak at (%d,%d)=%f' % (r0_last, z0_last, event_frame[r0_last, z0_last])

#        print 'Amplitude: %f Radial velocity: %f m/s' % ( amp.mean(), 0.01*(rzcom[0,-1] - rzcom[0,0]) / float(tau_a+tau_b)/dt )
#        # Plot poloidal vs. radial position
#        blob_t = np.arange(t0-tau_b, t0+tau_a)
#        plt.figure()
#        plt.subplot(121)
#        plt.plot(rzcom[0,:], rzcom[1,:])
#        plt.xlabel('R / cm')
#        plt.ylabel('z / cm')
#
#
#        # Plot blob size as a function of time
#        plt.subplot(122)
#        plt.plot(blob_t, fwhm_Rz[0,:], label='rad')
#        plt.plot(blob_t, fwhm_Rz[1,:], label='pol')
#        plt.plot(blob_t, amp, label='Intensity')
#        plt.xlabel('time / frames')
#        plt.ylabel('size / cm')
#        plt.legend()

        blob_amp[idx] = amp.mean()
        blob_vel[idx] = 0.01*(rzcom[0,-1] - rzcom[0,0]) / float(tau_a+tau_b)/dt
        blob_lp[idx]  = fwhm_Rz[0,:].mean()
        blob_lr[idx]  = fwhm_Rz[1,:].mean()

#    except:
#        event_ctr[idx] = -1
#        print 'Oooops'

plt.figure()
plt.plot(blob_amp, blob_vel, 'o')
plt.title('%d blob events' % (event_ctr == 1).sum() )
plt.xlabel('Amplitude / a.u.')
plt.ylabel('Velocity / ms^-1')

plt.figure()
plt.plot(blob_lp, blob_vel, 'o')
plt.title('%d blob events' % (event_ctr == 1).sum() )
plt.xlabel('Poloidal length / cm')
plt.ylabel('Velocity / ms^-1')

plt.figure()
plt.plot(blob_lr, blob_vel, 'o')
plt.title('%d blob events' % (event_ctr == 1).sum() )
plt.xlabel('Radial length / cm')
plt.ylabel('Velocity / ms^-1')


print 'Analyzed shot %d, averaged particle density %f' % ( shotnr, 0.1 )#shotparams['n_avg'] )
print '%d blobs analyzed, %d blobs ignored' % ( (event_ctr == 1).sum(), (event_ctr == -1).sum() )


plt.show()
