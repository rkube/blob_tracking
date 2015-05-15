#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

"""
=========
plotting
=========


Plotting routines for blobtrail objects


"""

import numpy as np
import matplotlib.pyplot as plt

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
        plt.xlabel('R / cm')
        plt.ylabel('Z / cm')

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



def plot_trail(blobtrail, frames, rz_array=None, xyi=None, trigger_box=None,
               sep_data=None, plot_com=False, plot_max=False,
               plot_shape=False, plot_geom=False, save_frames=False):
    """
    Plot the motion of the blob. The GPI frames are to be supplied
    externally

    Input:
        frames:         GPI data
        plot_com:       Mark the center of mass of the blob
        plot_max:       Mark the maximum of the blob
        plot_shape:     If available, mark the FWHM of the blob
        plot_geom:      Overplot triggering blox, limiter shadow and
                        separatrix
        save_frames:    Save the frames

    """

    # find the frames needed to track the object
    frame_idx = blobtrail.get_event_frames()
    print 'Plotting the blob event from frame %d-%d' % (frame_idx[0], frame_idx[-1])


    #    (self.event[1] + self.frame0 - self.tau_b,
    #     self.event[1] + self.frame0 + self.tau_f)
    #minval = (frames[self.event[1] + self.frame0 - self.tau_b:
    #                 self.event[1] + self.frame0 + self.tau_f, :, :]).min()
    #maxval = (frames[self.event[1] + self.frame0 - self.tau_b:
    #                 self.event[1] + self.frame0 + self.tau_f, :, :]).max()
    #frames[self.event[1] + self.frame0 - self.tau_b:
    #       self.event[1] + self.frame0 + self.tau_f, 0, 0] = minval
    #frames[self.event[1] + self.frame0 - self.tau_b:
    #       self.event[1] + self.frame0 + self.tau_f, 0, 1] = maxval

    # Get minimum and maximum value for color scales
    minval = frames[frame_idx, :, :].min()
    maxval = frames[frame_idx, :, :].max()
    print 'min = %f, max = %f' % (minval, maxval)
    print 'plotting from %d - %d' % (self.tau_b, self.tau_f)

    num_levels = 64
    color_levels = np.linspace(minval, maxval, num_levels)

    for f_idx, tau in enumerate(blobtrail.get_tau()):
        plt.figure()
        plt.title('frame %05d' % (self.event[1] + self.frame0 + tau))
        plt.xlabel('R / cm')
        plt.ylabel('Z / cm')

        # Try plotting everythin in machine coordinates. If it fails,
        # draw in pixels
        try:
            1 / 0
            zi = griddata(rz_array.reshape(64*64, 2),
                          frames[self.event[1] +
                                 self.frame0 + tau, :, :].reshape(64*64),
                          xyi.reshape(64*64, 2), method='linear')
#                zi[0] = 5.0#np.max(frames)
#                zi[1] = 5.0#np.max(frames)
            plt.contour(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64),
                        32, linewidths=0.5, colors='k')
            plt.contourf(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64),
                         32, cmap=plt.cm.hot, levels=np.linspace(0.0,
                                                                 maxval,
                                                                 32))

        except:
            plt.contour(frames[self.event[1] + self.frame0 + tau, :, :],
                        22, linewidths=0.5, colors='k')
            plt.contourf(frames[self.event[1] + self.frame0 + tau, :, :],
                         32, cmap=plt.cm.hot, levels=np.linspace(0.0,
                                                                 maxval,
                                                                 32))

        plt.colorbar(ticks=np.arange(0.0, 3.5, 0.5), format='%3.1f')

        if plot_com:
            try:
                1/0
                if plot_shape:
                    frame_xerr = self.fwhm_ell_rad[:f_idx+1]
                    frame_xerr[:-1] = 0.
                    frame_yerr = self.fwhm_ell_pol[:f_idx+1]
                    frame_yerr[:-1] = 0.
                    plt.errorbar(xyi[self.xycom[:f_idx + 1,
                                                0].astype('int'),
                                     self.xycom[:f_idx + 1,
                                                1].astype('int'), 0],
                                 xyi[self.xycom[:f_idx + 1,
                                                0].astype('int'),
                                     self.xycom[:f_idx + 1,
                                                1].astype('int'), 1],
                                 xerr=frame_xerr, yerr=frame_yerr,
                                 ecolor='w', linestyle='None',
                                 mfc='white', mec='green', marker='s')

                else:
                    plt.plot(xyi[self.xycom[:f_idx + 1, 0].astype('int'),
                                 self.xycom[:f_idx + 1, 1].astype('int'),
                                 0],
                             xyi[self.xycom[:f_idx + 1, 0].astype('int'),
                                 self.xycom[:f_idx + 1, 1].astype('int'),
                                 1], '-ws')

                # Set the coordinates for plotting the text field
                text_x, text_y = 86.2, -6.
            except:
                plt.plot(self.xycom[:f_idx + 1, 1],
                         self.xycom[:f_idx + 1, 0], '-bs')
                text_x, text_y = 5., 2.

            #if (tau < self.tau_f - 1):
            #    plt.text(text_x, text_y, '$V_{COM} = (%4.1f, %4.1f)$' %
            #             (self.get_velocity_com(rz_array)[f_idx, 0],
            #              self.get_velocity_com(rz_array)[f_idx, 1]),
            #             fontdict=dict(size=16., color='white',
            #                           weight='bold'))

        if plot_max:
            try:
                plt.plot(xyi[self.xymax[:f_idx + 1, 0],
                             self.xymax[:f_idx + 1, 1], 0],
                         xyi[self.xymax[:f_idx + 1, 0],
                             self.xymax[:f_idx + 1, 1], 1], '-.wo')
                text_x, text_y = 86.2, -6.

            except TypeError:
                plt.plot(self.xymax[:f_idx+1, 1], self.xymax[:f_idx+1, 0],
                         '-.wo')
                text_x, text_y = 5., 2.

            if (tau < self.tau_f - 1):
                plt.text(text_x, text_y, '$V_{max} = (%4.1f, %4.1f)$' %
                         (self.get_velocity_max(rz_array)[f_idx, 0],
                          self.get_velocity_max(rz_array)[f_idx, 1]),
                         fontdict=dict(size=16.,
                                       color='white', weight='bold'))

        #if plot_geom:
        #    try:
        #        # Get the position of the pixels for the separatrix
        #        # and limiter
        #        separatrix_pxs = surface_line(sep_data['rmid'].
        #                                      reshape(64, 64) >
        #                                      sep_data['rmid_sepx'],
        #                                      mode='max')
        #        limiter_pxs = surface_line(sep_data['rmid'].
        #                                   reshape(64, 64) <
        #                                   sep_data['rmid_lim'],
        #                                   mode='min')

        #        # Compute position, width and height of the triggering box
        #        # tb_lower_left = (xyi[trigger_box[2], trigger_box[0], 0],
        #        #                  xyi[trigger_box[2], trigger_box[0], 1])
        #        # tb_width = (xyi[trigger_box[2], trigger_box[1], 0] -
        #        #             xyi[trigger_box[2], trigger_box[0], 0])
        #        # tb_height = (xyi[trigger_box[3], trigger_box[0], 1] -
        #        #              xyi[trigger_box[2], trigger_box[0], 1])

        #        # Plot the triggering domain. Position, height and width
        #        # are not automatically determined but static values.

        #        triggering_box = mpatches.Rectangle((89.9, -4.5),
        #                                            width=1.0, height=3.2,
        #                                            fill=False,
        #                                            ls='dashdot',
        #                                            ec='w', lw=3)
        #        fig = plt.gcf()
        #        ax = fig.gca()
        #        ax.add_patch(triggering_box)

        #        # Plot the separatrix
        #        sep_x = [xyi[i, separatrix_pxs[i], 0] for i in
        #                 np.arange(64)]
        #        sep_y = [xyi[i, separatrix_pxs[i], 1] for i in
        #                 np.arange(64)]
        #        plt.plot(sep_x, sep_y, 'w--', linewidth=4)

        #        lim_x = [xyi[i, limiter_pxs[i], 0] for i in np.arange(64)]
        #        lim_y = [xyi[i, limiter_pxs[i], 1] for i in np.arange(64)]
        #        plt.plot(lim_x, lim_y, 'w-.', linewidth=4)

        #    except:
        #        print 'Error plotting geometry :('

        if save_frames:
            F = plt.gcf()
            F.savefig('%d/frames/frame_%05d.eps' % (self.shotnr,
                                                    self.event[1] +
                                                    self.frame0 + tau))
            plt.close()

    plt.show()


# End of file plotting.py
