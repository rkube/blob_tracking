#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

"""
=========
plotting
=========


Plotting routines for blobtrail objects

* plot_trail_simple -> Simple plotting without special geometry
* plot_trail_geom -> Use outboard midplane geometry of CMOD. Needs geometry and separatrix data

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d
import matplotlib.patches as mpatches
from geometry import velocity_max, velocity_com

def plot_trail_simple(blobtrail, frames, plot_com='o', plot_max='^', plot_shape=True, save_frames=False):
    """
    Plot the motion of the blob. The GPI frames are to 
    be supplied externally

    Input:
        frames:         GPI data
        plot_com:       Mark the center of mass of the blob
        plot_max:       Mark the maximum of the blob
        plot_shape:     If available, plot the radial extend of the blob
        save_frames:    Save the frames

    """

    # find the frames needed to track the object
    frame_idx = blobtrail.get_event_frames()

    # Length of the blobtrail object
    tau = np.arange(blobtrail.get_tau().size)

    # Get the coordinates of the blob
    xymax = blobtrail.get_xymax()
    xycom = blobtrail.get_xycom()

    ell_rad = blobtrail.get_ell_rad()
    ell_pol = blobtrail.get_ell_pol()

    # Get minimum and maximum value for color scales
    minval = frames[frame_idx, :, :].min()
    maxval = frames[frame_idx, :, :].max()

    num_levels = 64
    color_levels = np.linspace(minval, maxval, num_levels)

    for fidx, tau_idx in zip(frame_idx, tau):
        plt.figure()
        plt.title('frame %05d' % (fidx))
        plt.xlabel('x / px')
        plt.ylabel('y / px')

        plt.contourf(frames[fidx], cmap=plt.cm.hot, levels=color_levels)

        if plot_com is not None:
            plt.plot(xycom[:tau_idx, 1], xycom[:tau_idx, 0], plot_com)

        # Plot maximum with errorbars
        if plot_max is not None:
            if plot_shape:
                plt.errorbar(xymax[:tau_idx, 1], xymax[:tau_idx, 0], 
                             xerr=ell_rad[:tau_idx], yerr=ell_pol[:tau_idx],
                             ecolor='w', linestyle='None',
                             mfc='white', mec='green', marker=plot_max)
            else:
                plt.plot(xymax[:tau_idx, 1], xymax[:tau_idx, 0], plot_max)

        plt.colorbar()

    if save_frames:
        F = plt.gcf()
        F.savefig('%d/frames/frame_%05d.eps' % (self.shotnr,
                                                self.event[1] +
                                                self.frame0 + tau))
        plt.close()

    plt.show()



def plot_trail_geom(blobtrail, frames, rz_array=None, xyi=None, trigger_box=None,
                    sep_data=None, plot_com=False, plot_max=False,
                    plot_shape=False, plot_geom=False, save_frames=False):
    """
    Plot the motion of the blob using CMOD GPI geometry
    The GPI frames are to be supplied externally.

    Input:
        frames:         ndarray, GPI data. axis0: time, axis1: poloidal, axis2: radial
        rz_array:       ndarray, Where the GPI data is known
        xyi:            ndarray, Array on which we interpolate GPI data for output
        trigger_box:    ndarray, Array where blobs have been detected
        sep_data:       IDL .sav structure, see /home/rkube/IDL/separatrix.pro
        plot_com:       Mark the center of mass of the blob
        plot_max:       Mark the maximum of the blob
        plot_shape:     If available, mark the FWHM of the blob
        plot_geom:      Overplot triggering blox, limiter shadow and
                        separatrix
        save_frames:    Save the frames
    """

    # find the frames needed to track the object
    frame_idx = blobtrail.get_event_frames()
    tau = np.arange(blobtrail.get_tau().size)
    xymax = blobtrail.get_xymax()
    xycom = blobtrail.get_xycom()
    ell_rad = blobtrail.get_ell_rad()
    ell_pol = blobtrail.get_ell_pol()

    # Get minimum and maximum value for color scales
    minval = frames[frame_idx, :, :].min()
    maxval = frames[frame_idx, :, :].max()
    print 'min = %f, max = %f' % (minval, maxval)
    print 'plotting from %d - %d' % (tau[0], tau[-1])

    num_levels = 64
    color_levels = np.linspace(minval, maxval, num_levels)

    for fidx, tau_idx in zip(frame_idx, tau):
        plt.figure()
        plt.title('frame %05d' % (fidx))
        plt.xlabel('R / cm')
        plt.ylabel('Z / cm')

        # Try plotting everythin in machine coordinates. If it fails draw in pixels
        #try:
        zi = griddata(rz_array.reshape(64 * 64, 2),
                frames[fidx, :, :].reshape(64 * 64),
                xyi.reshape(64 * 64, 2), method='linear')
        #plt.contour(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64),
        #            32, linewidths=0.5, colors='k')
        plt.contourf(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64),
                     num_levels, cmap=plt.cm.hot, levels=color_levels)

        #except:
        #    print 'Failed to grid data on rz_array... Drawing in pixels'
        #    #plt.contour(frames[self.event[1] + self.frame0 + tau, :, :],
        #    #            32, linewidths=0.5, colors='k')
        #    plt.contourf(frames[fidx, :, :], num_levels, cmap=plt.cm.hot, levels=color_levels)
        plt.colorbar(ticks=np.arange(minval, maxval, (maxval - minval) / 5.), format='%3.1f')

        if plot_com:
            # Plot the leading up blob trail without error bars
            plt.plot(xyi[xycom[:tau_idx + 1, 0].astype('int'),
                         xycom[:tau_idx + 1, 1].astype('int'),
                         0],
                     xyi[xycom[:tau_idx + 1, 0].astype('int'),
                         xycom[:tau_idx + 1, 1].astype('int'),
                     1], '-ws')

            # Interpolate width from pixel to physical units
            ip_rad = interp1d(np.arange(64), xyi[xycom[tau_idx, 0].astype('int'), :, 0], kind='quadratic')
            ip_pol = interp1d(np.arange(64), xyi[:, xycom[tau_idx, 1].astype('int'), 1], kind='quadratic')
            frame_xerr = ip_rad(np.array([xycom[tau_idx, 0] - ell_rad[tau_idx], xycom[tau_idx, 0] + ell_rad[tau_idx]]))
            frame_yerr = ip_pol(np.array([xycom[tau_idx, 0] - ell_pol[tau_idx], xycom[tau_idx, 0] + ell_pol[tau_idx]]))

            frame_xerr = np.abs(frame_xerr[1] - frame_xerr[0])
            frame_yerr = np.abs(frame_yerr[1] - frame_yerr[0])
            # Plot current blob position with error bars
            plt.errorbar(xyi[xycom[tau_idx, 0].astype('int'),
                             xycom[tau_idx, 1].astype('int'), 0],
                         xyi[xycom[tau_idx, 0].astype('int'),
                             xycom[tau_idx, 1].astype('int'), 1],
                         xerr=frame_xerr, yerr=frame_yerr, ecolor='w', linestyle='None', mfc='white', mec='green', marker='s')

            # Set the coordinates for plotting the text field
            text_x, text_y = 86.2, -6.

            if (tau_idx > 0):
                plt.text(text_x, text_y, '$V_{COM} = (%4.1f, %4.1f)$' %
                         (velocity_com(blobtrail, rz_array)[tau_idx - 1, 0],
                          velocity_com(blobtrail, rz_array)[tau_idx - 1, 1]),
                          fontdict=dict(size=16., color='white',
                                        weight='bold'))

        if plot_max:
            plt.plot(xyi[xymax[:tau_idx + 1, 0].astype('int'),
                         xymax[:tau_idx + 1, 1].astype('int'), 0],
                     xyi[xymax[:tau_idx + 1, 0].astype('int'),
                         xymax[:tau_idx + 1, 1].astype('int'), 1], '-.wo')
            text_x, text_y = 86.2, -6.

            if (tau_idx > 0): 
                vmax_str = r"$V_{max} = (%4.1f, %4.1f)$" % (velocity_max(blobtrail, rz_array=rz_array)[tau_idx - 1, 0],
                                                            velocity_max(blobtrail, rz_array=rz_array)[tau_idx - 1, 1]),
                plt.text(text_x, text_y, vmax_str, 
                         fontdict=dict(size=16., color='white', weight='bold'))

        if plot_geom:
            # Get the position of the pixels for the separatrix
            # and limiter
            separatrix_pxs = surface_line(sep_data['rmid'].
                                          reshape(64, 64) >
                                          sep_data['rmid_sepx'],
                                          mode='max')
            limiter_pxs = surface_line(sep_data['rmid'].
                                       reshape(64, 64) <
                                       sep_data['rmid_lim'],
                                       mode='min')

                # Compute position, width and height of the triggering box
                # tb_lower_left = (xyi[trigger_box[2], trigger_box[0], 0],
                #                  xyi[trigger_box[2], trigger_box[0], 1])
                # tb_width = (xyi[trigger_box[2], trigger_box[1], 0] -
                #             xyi[trigger_box[2], trigger_box[0], 0])
                # tb_height = (xyi[trigger_box[3], trigger_box[0], 1] -
                #              xyi[trigger_box[2], trigger_box[0], 1])

                # Plot the triggering domain. Position, height and width
                # are not automatically determined but static values.

            triggering_box = mpatches.Rectangle((89.9, -4.5),
                                                width=1.0, height=3.2,
                                                fill=False,
                                                ls='dashdot',
                                                ec='w', lw=3)
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_patch(triggering_box)

            # Plot the separatrix
            sep_x = [xyi[i, separatrix_pxs[i], 0] for i in
                     np.arange(64)]
            sep_y = [xyi[i, separatrix_pxs[i], 1] for i in
                     np.arange(64)]
            plt.plot(sep_x, sep_y, 'w--', linewidth=4)

            lim_x = [xyi[i, limiter_pxs[i], 0] for i in np.arange(64)]
            lim_y = [xyi[i, limiter_pxs[i], 1] for i in np.arange(64)]
            plt.plot(lim_x, lim_y, 'w-.', linewidth=4)


        if save_frames:
            F = plt.gcf()
            F.savefig('%d/frames/frame_%05d.eps' % (self.shotnr,
                                                    self.event[1] +
                                                    self.frame0 + tau))
            plt.close()

    plt.show()


def surface_line(sep_pixels, mode='max'):
    """
    Given the pixels which are mapped to the closed field line region,
    return the pixel with the largest radial coordinate for each
    poloidal coordinate.
    """

    # Index all pixels radially
    lin_array = np.repeat(np.arange(64), 64).reshape(64, 64).T

    # Apply the mask
    la_masked = np.ma.array(lin_array, mask=sep_pixels)

    # Return the maximum radial indices for each poloidal position
    if (mode == 'max'):
        return la_masked.argmax(axis=1)
    elif (mode == 'min'):
        return la_masked.argmin(axis=1)


# End of file plotting.py
