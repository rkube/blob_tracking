#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import blobtrail
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import helper_functions
import geometry 


def velocity_analysis(trails, frames, sol_px, rz_array, xyi):
    """
    Study blob velocity dependence on cross-field size

    Input:
        trails: list, List of blob trail events
        frames: ndarray, axis0: time
        sol_px: ndarray, index of pixels that are in the SOL
    """
    
    # List of pixels that are in the SOL
    sol_px_list = sol_px.tolist()

    # Average cross-field size of the blob while in the SOL
    blob_ell_rad = np.zeros([len(trails)])
    blob_ell_pol = np.zeros([len(trails)])

    blob_vcom_rad = np.zeros([len(trails)])
    blob_vcom_pol = np.zeros([len(trails)])

    # Number of blobs we have analyzed
    blob_count = 0
    for idx, trail in enumerate(trails):
        print 'trail %d / %d' % (idx, len(trails))
        # Find the instances, where the current blobtrail is recorded
        # in the scrape-off layer
        good_pos_idx = geometry.blob_in_sol(trail, sol_px_list, logger=None)
        if ( good_pos_idx.sum() < 5 ):
            continue
        blob_count += 1

#        plt.figure()
#        plt.contourf(frames[trail.get_event()[1], :, :])
#        plt.plot(sol_px[:, 1], sol_px[:, 0], 'k.')
#        plt.plot(trail.get_xycom().astype('int')[:, 1], trail.get_xycom().astype('int')[:, 0], 'ro')
#        plt.show()

        # Determine mean blob size in SOL
        xycom = trail.get_xycom()
        ell_rad_px = trail.get_ell_rad()
        ell_pol_px = trail.get_ell_pol()
        ell_rad = np.zeros_like(ell_rad_px)
        ell_pol = np.zeros_like(ell_pol_px)

        # Interpolate the width, given by ell_rad and ell_pol on the physical grid
        for tau_idx, tau in enumerate(trail.get_tau()):
            ip_rad = interp1d(np.arange(64), xyi[xycom[tau_idx, 0].astype('int'), :, 0], kind='quadratic')
            ip_pol = interp1d(np.arange(64), xyi[:, xycom[tau_idx, 1].astype('int'), 1], kind='quadratic')
            try:
                tau_xerr = ip_rad(np.array([xycom[tau_idx, 0] - ell_rad_px[tau_idx], xycom[tau_idx, 0] + ell_rad_px[tau_idx]]))
                ell_rad[tau_idx] = np.abs(tau_xerr[1] - tau_xerr[0])
            except ValueError:
                ell_rad[tau_idx] = ell_rad[tau_idx - 1]

            try:
                tau_yerr = ip_pol(np.array([xycom[tau_idx, 0] - ell_pol_px[tau_idx], xycom[tau_idx, 0] + ell_pol_px[tau_idx]]))
                ell_pol[tau_idx] = np.abs(tau_yerr[1] - tau_yerr[0])
            except ValueError:
                ell_pol[tau_idx] = ell_pol[tau_idx - 1]

        blob_ell_rad[idx] = ell_rad[good_pos_idx].mean()
        blob_ell_pol[idx] = ell_pol[good_pos_idx].mean()

        # Compute average blob velocity
        # We compute the blob velocity with a centered difference scheme.
        # Thus, when indexing the velocity with good_pos_idx, we have to discard
        # the first and last position
        vcom = geometry.velocity_com(trail, rz_array)
        print 'mean(Vcom):rad=%f, pol=%f' % (vcom.mean(axis=0)[0], vcom.mean(axis=0)[1])

        blob_vcom_rad[idx] = vcom[good_pos_idx[1:]].mean(axis=0)[0]
        blob_vcom_pol[idx] = vcom[good_pos_idx[1:]].mean(axis=0)[1]



    title_str = "%d trails" % (len(trails))

    fig = plt.figure(figsize=(8, 12))
    fig.text(0.5, 0.95, title_str, ha='center')
    plt.subplot(411)
    plt.hist(blob_ell_rad)
    plt.ylabel(r"$\ell_{\mathrm{rad}} / \mathrm{cm}$")

    plt.subplot(412)
    plt.hist(blob_ell_pol)
    plt.ylabel(r"$\ell_{\mathrm{pol}} / \mathrm{cm}$")

    plt.subplot(413)
    plt.hist(blob_vcom_rad)
    plt.title('blob_vcom_rad')
    plt.ylabel(r"$V_{\mathrm{rad}} / \mathrm{ms}^{-1}$")

    plt.subplot(414)
    plt.hist(blob_vcom_pol)
    plt.ylabel(r"$V_{\mathrm{pol}} / \mathrm{ms}^{-1}$")

    plt.show()

#
#        # Interpolate velocity on x_sol
#        f = interp1d( trail.get_trail_com()[good_pos_idx, 1], trail.get_velocity_com()[good_pos_idx, 1], bounds_error = False)
#        vel_ip = f(x_sol)
#
#        if ( blob_ell[idx,0] < 0.4 ):
#            count_idx_low += np.invert( np.isnan(vel_ip) )
#            vel_ip[ np.isnan(vel_ip) ] = 0.0
#            mean_v_low += vel_ip        
#
#        elif ( blob_ell[idx,0] > 0.4 and blob_ell[idx,0] < 0.6 ):
#            count_idx_med += np.invert( np.isnan(vel_ip) )
#            vel_ip[ np.isnan(vel_ip) ] = 0.0
#            mean_v_med += vel_ip        
#
#        else:
#            count_idx_large += np.invert( np.isnan(vel_ip) )
#            vel_ip[ np.isnan(vel_ip) ] = 0.0
#            mean_v_large += vel_ip        
#
#    print count_idx_low
#    print count_idx_med
#    print count_idx_large
#
#    print 'Accepted %d blobs' % ( blob_count )
#
#    plt.figure()
#    plt.plot(x_sol, mean_v_low / count_idx_low, label='Low' )
#    plt.plot(x_sol, mean_v_med / count_idx_med, label='Medium' )
#    plt.plot(x_sol, mean_v_large / count_idx_large, label='Large' )
#    plt.legend()
#    plt.show()
#
#

# End of file velocity_analysis.py
