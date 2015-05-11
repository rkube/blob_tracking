#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-


import numpy as np
import matplotlib as mpl
# mpl.use('AGG')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os.path
from scipy.io import readsav
from misc.load_mdsframes import load_mdsframes
from misc.phantom_helper import make_rz_array
from scipy.interpolate import griddata
from plotting.separatrix_line import surface_line

shotnr = 1100803011
# frame0 = 20000
frames, fi = load_mdsframes(shotnr, test=False)
f_max = np.max(frames)
f_min = np.min(frames)
# frames = frames[frame0:25000,:,:]

# Compute the indices between LCFS and limiter shadow
s = readsav('%d/%d_separatrix.sav' % (shotnr, shotnr), verbose=False)
separatrix_pxs = surface_line(s['rmid'].reshape(64, 64) > s['rmid_sepx'],
                              mode='max')
limiter_pxs = surface_line(s['rmid'].reshape(64, 64) < s['rmid_lim'],
                           mode='min')

rz_array, transform_data = make_rz_array(fi)
xxi, yyi = np.meshgrid(np.linspace(rz_array[:, :, 0].min(),
                                   rz_array[:, :, 0].max(), 64),
                       np.linspace(rz_array[:, :, 1].min(),
                                   rz_array[:, :, 1].max(), 64))
xyi = np.concatenate((xxi[:, :, np.newaxis], yyi[:, :, np.newaxis]), axis=2)
maxxxx = np.max(frames)


for frame in np.arange(17000, 20000):
    if os.path.isfile('%d/frames/frame_%05d.png' % (shotnr, frame)):
        print 'Skipping frame %d' % (frame)
        continue
    print 'Plotting frame %d' % (frame)

    fig = plt.figure()
    zi = griddata(rz_array.reshape(64 * 64, 2),
                  frames[frame, :, :].reshape(64 * 64),
                  xyi.reshape(64 * 64, 2), method='linear')
    zi[0] = 5.0
    zi[1] = 0.5
    plt.contour(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64), 32,
                linewidths=0.5, colors='k')
    plt.contourf(xyi[:, :, 0], xyi[:, :, 1], zi.reshape(64, 64), 32,
                 cmap=plt.cm.hot)
    plt.colorbar(ticks=np.arange(0.0, maxxxx, 0.5), format='%3.1f')
    plt.title('frame %d' % frame)
    plt.xlabel('R / cm')
    plt.ylabel('Z / cm')

    print 'Trigger box: lower left: R=%f Z=%f, width = %f, height = %f' % \
        (xyi[16, 40, 0], xyi[16, 40, 1],
         (xyi[16, 50, 0] - xyi[16, 40, 0]),
         (xyi[48, 40, 1] - xyi[16, 40, 1]))
    tb_lower_left = (xyi[16, 40, 0], xyi[16, 40, 1])
    tb_width = (xyi[16, 50, 0] - xyi[16, 40, 0])
    tb_height = (xyi[48, 40, 1] - xyi[16, 40, 1])

    # Plot the triggering domain
    triggering_box = mpatches.Rectangle(tb_lower_left, width=tb_width,
                                        height=tb_height, fill=False,
                                        ls='dashdot',
                                        ec='w', lw=3)
    ax = fig.gca()
    ax.add_patch(triggering_box)

    # Plot the separatrix
    sep_x = [xyi[i, separatrix_pxs[i], 0] for i in np.arange(64)]
    sep_y = [xyi[i, separatrix_pxs[i], 1] for i in np.arange(64)]
    plt.plot(sep_x, sep_y, 'w--', linewidth=4)

    lim_x = [xyi[i, limiter_pxs[i], 0] for i in np.arange(64)]
    lim_y = [xyi[i, limiter_pxs[i], 1] for i in np.arange(64)]
    plt.plot(lim_x, lim_y, 'w-.', linewidth=4)

    F = plt.gcf()
    F.savefig('%d/frames/frame_%05d.png' % (shotnr, frame))
    plt.close()

plt.show()

# End of file plot_movie_frames.py
