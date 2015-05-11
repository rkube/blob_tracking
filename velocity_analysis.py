#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import blobtrail
import cPickle
import helper_functions
from load_mdsframes import load_mdsframes

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

"""
Study the velocity dependence on blob length
"""

shotnr = 1100803020
thresh = 2.5


picklefile = open('%d/%d_trails_%2d.pkl' % (shotnr, shotnr, int(10*thresh) ), 'rb' )
trails = cPickle.load(picklefile)
picklefile.close()
frames, fi = load_mdsframes(shotnr, test = True )
frames = 0.0
    
# Find the SOL and the corresponding indices
mask = helper_functions.find_sol_mask(shotnr)
good_idx = helper_functions.find_sol_pixels(shotnr, fi)
domain = np.array(good_idx)[:,:].tolist()


blob_ell = np.zeros([len(trails), 2])
blob_count = 0


x_sol = np.arange(36,50,0.1)
mean_v_low = np.zeros_like(x_sol)
count_idx_low = np.zeros_like(x_sol)

mean_v_med = np.zeros_like(x_sol)
count_idx_med = np.zeros_like(x_sol)

mean_v_large = np.zeros_like(x_sol)
count_idx_large = np.zeros_like(x_sol)


for idx, trail in enumerate(trails[:]):
    good_pos_idx = helper_functions.blob_in_sol(trail, domain, logger = None)
    if ( np.sum(good_pos_idx) < 10 ):
        continue
    blob_count += 1

    # Determine mean blob size in SOL
    blob_ell[idx,0] = trail.get_ell_rad()[good_pos_idx].mean()#[[newtrail.get_frame0()]]

    # Interpolate velocity on x_sol
    f = interp1d( trail.get_trail_com()[good_pos_idx, 1], trail.get_velocity_com()[good_pos_idx, 1], bounds_error = False)
    vel_ip = f(x_sol)

    if ( blob_ell[idx,0] < 0.4 ):
        count_idx_low += np.invert( np.isnan(vel_ip) )
        vel_ip[ np.isnan(vel_ip) ] = 0.0
        mean_v_low += vel_ip        

    elif ( blob_ell[idx,0] > 0.4 and blob_ell[idx,0] < 0.6 ):
        count_idx_med += np.invert( np.isnan(vel_ip) )
        vel_ip[ np.isnan(vel_ip) ] = 0.0
        mean_v_med += vel_ip        

    else:
        count_idx_large += np.invert( np.isnan(vel_ip) )
        vel_ip[ np.isnan(vel_ip) ] = 0.0
        mean_v_large += vel_ip        

print count_idx_low
print count_idx_med
print count_idx_large

print 'Accepted %d blobs' % ( blob_count )

plt.figure()
plt.plot(x_sol, mean_v_low / count_idx_low, label='Low' )
plt.plot(x_sol, mean_v_med / count_idx_med, label='Medium' )
plt.plot(x_sol, mean_v_large / count_idx_large, label='Large' )
plt.legend()
plt.show()
