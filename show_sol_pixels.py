#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from misc import helper_functions, load_mdsframes
from misc.phantom_helper import make_rz_array


"""
Mark the pixels within the SOL for a given shot
"""


shotnr = 1120711005
frames, frames_info, frames_mean, frames_rms = load_mdsframes.load_mdsframes_mean_rms( shotnr, test = True )

good_idx = helper_functions.find_sol_pixels(shotnr, frames_info)
domain = np.array(good_idx)[:,:].tolist()

print frames_info

plt.figure()
plt.title('# %d' % shotnr )

rz_array, transform_data = make_rz_array(frames_info)
xxi, yyi = np.meshgrid( np.linspace( np.min(rz_array[:,:,0] ), np.max( rz_array[:,:,0] ),64 ), np.linspace( np.min( rz_array[:,:,1] ), np.max( rz_array[:,:,1] ),64 ) )
xyi = np.concatenate( (xxi[:,:,np.newaxis], yyi[:,:,np.newaxis]), axis=2 )

for px in domain:
    plt.plot( xyi[px[0], px[1], 0], xyi[px[0], px[1], 1], 'k.' )

plt.contour(xyi[:,:,0], xyi[:,:,1], frames_mean, 32, linewidths = 0.5, colors = 'k')
plt.contourf(xyi[:,:,0], xyi[:,:,1], frames_mean, 32, cmap = plt.cm.hot)
plt.colorbar(ticks=np.arange(0.0, np.max(frames_mean), 25.0), format='%3.1f')

plt.show()
