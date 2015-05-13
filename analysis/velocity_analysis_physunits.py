#!/opt/local/bin/python
#-*- encoding: utf-8 -*-

import numpy as np
from plotting.figure_defs import set_rcparams_pres_small
import matplotlib as mpl
set_rcparams_pres_small(mpl.rcParams)
import matplotlib.pyplot as plt
from blob_tracking.blobtrail import blobtrail
import cPickle
from misc import helper_functions
from misc.phantom_helper import make_rz_array
from misc.load_mdsframes import load_mdsframes


"""
Plot blob velocities for each shot in physical units.
Superpose sheath velocity scaling
"""


# Define physical parameters for this show

R = 0.68                                # Major radius in m
a = 0.22                                # Minor radius in m
q = 1.602e-19                           # Elementary charge
m_p = 1.67e-27                          # Proton mass
k_b = 1.38e-23                          # Boltzmann constant
L_conn = 5.
Z = 2                                   # We have deuterium ions
T_e = 10.                               # Electron temperature in eV
T_i = 10.                               # Ion temperature in eV
m_i = Z * m_p                           # Ion mass ( deuterium )
B = 4.0 / ( 1. + a/R )                  # Magnetic field in SOL, in T
omega_ci =  q*B/m_i                     # Ion cyclotron frequency
cs = np.sqrt(q * (T_e+T_i) / m_i )      # Acoustic velocity
rho_s = cs / omega_ci                   # Ion thermal gyroradius
ell_star = (2. * rho_s**4. * L_conn**2 / R )**(1./5.)
v_star = cs * np.sqrt(2. * ell_star / R)

print 'v_star = %f m/s' % ( v_star )
print 'ell_star = %f m' % ( ell_star )


# Analytic blob velocity scaling
v_analytic = lambda l, dn, c1, c2: c2/2. * ell**3 * ( -1. + np.sqrt(1. + 4/ell**5 * c1/c2**2 * (dn/(1.+dn)) ) )
v_small = lambda l, dn, c1, c2: np.sqrt(l*(dn/(1.+dn))*c1)
datafile = np.load('/Users/ralph/Dropbox/Kube/blob-sheath/graphics/data/sheath_c12.npz')
print datafile['dtlist']
c1_list = datafile['c1']
c2_list = datafile['c2']

shotnr = 1120217008
n = 0.15

# Load blob trails for the shot
thresh = 1.5
#picklefile = open('1120217021/1120217021_trails_thresh15_20120410.pkl', 'rb' )
picklefile = open('scan_120217/trails_ne%2d_thresh15_gaussian_20120424.pkl' % ( int(100*n) ), 'rb')
trails = cPickle.load(picklefile)
picklefile.close()
frames, frame_info = load_mdsframes(shotnr, test = True )
frames = 0.0

# Find the SOL and the corresponding indices
mask = helper_functions.find_sol_mask(shotnr)
good_idx = helper_functions.find_sol_pixels(shotnr, frame_info)
domain = np.array(good_idx)[:,:].tolist()


# Use physical units
rz_array, transform_data = make_rz_array(frame_info)

#blob_ell_s = np.zeros([len(trails)])
#blob_vel_s = np.zeros([len(trails)])
#
#blob_ell_m = np.zeros([len(trails)])
#blob_vel_m = np.zeros([len(trails)])
#
#blob_ell_l = np.zeros([len(trails)])
#blob_vel_l = np.zeros([len(trails)])

blob_ell = np.zeros(len(trails))
blob_vel = np.zeros(len(trails))

blob_count = 0

# Plot analytic blob velocity scaling
ell = np.arange(0,5,0.01)
fig = plt.figure( figsize = (6,5))
fig.subplots_adjust(bottom = 0.15, left = 0.15)

# If unicode plotting ever works on a mac
#plt.plot(ell, v_analytic(ell, 0.1,   c1_list[0], c2_list[0]), 'k', linewidth=3., label=u'∆n/N = 0.1')
#plt.plot(ell, v_analytic(ell, 1.0,   c1_list[4], c2_list[4]), 'k-.', linewidth=3., label=u'∆n/N = 1.0')
#plt.plot(ell, v_analytic(ell, 10.0,  c1_list[6], c2_list[6]), 'k-.', linewidth=3., label=u'∆n/N = 5.0')
#plt.xlabel(u'ℓ/ℓ_*', fontproperties = font)
#plt.ylabel(u'V/V_*', fontproperties = font)

v_analytic_small =   v_analytic(ell, 0.1,   c1_list[0], c2_list[0])
v_analytic_medium =  v_analytic(ell, 1.0,   c1_list[4], c2_list[4])
v_analytic_large =   v_analytic(ell, 10.0,  c1_list[6], c2_list[6])

plt.plot(ell, v_analytic_small, 'k-.', linewidth=1., label='$\\mathbf{\\Delta n/N = 0.1}$')
plt.plot(ell, v_analytic_medium, 'k', linewidth=3., label='$\\mathbf{\\Delta n/N = 1.0}$')
plt.plot(ell, v_analytic_large, 'k-.', linewidth=1., label='$\\mathbf{\\Delta n/N = 5.0}$')
#plt.plot(ell, v_small   (ell, 0.1,   c1_list[0], c2_list[0]), 'r-.', linewidth=3., label='')
plt.plot(ell, v_small   (ell, 1.0,   c1_list[4], c2_list[4]), 'r-', linewidth=3., label='')
#plt.plot(ell, v_small   (ell, 10.0,  c1_list[6], c2_list[6]), 'r--', linewidth=3., label='')
plt.xlabel('$\\ell / \\ell_*$', fontsize=30)
plt.ylabel('$V / V_*$', fontsize=30)


ax = plt.gca()
ax.fill_between(ell, v_analytic_small, v_analytic_large, facecolor='blue', alpha=0.5)



for idx, trail in enumerate(trails):
    good_pos_idx = helper_functions.blob_in_sol(trail, domain, logger = None)
    
    if ( np.sum(good_pos_idx) < 8 ):
        continue
    blob_count += 1
    
    # Determine mean blob size in SOL
    size = trail.get_ell_rad()[good_pos_idx]

    blobamp = trail.get_amp().mean()

    blob_ell[idx] = size[size>0].mean()*1e-2
    blob_vel[idx] = trail.get_velocity_com(rz_array)[good_pos_idx, 0].mean()

#    if ( blobamp < 3.0 ):
#        blob_ell_s[idx] = size[size>0].mean()*1e-2
#        # Determine mean radial blob velocity in SOL
#        blob_vel_s[idx] = trail.get_velocity_com(rz_array)[good_pos_idx, 0].mean()
#
#    elif ( blobamp >= 3.0 ):
#        blob_ell_m[idx] = size[size>0].mean()*1e-2
#        # Determine mean radial blob velocity in SOL
#        blob_vel_m[idx] = trail.get_velocity_com(rz_array)[good_pos_idx, 0].mean()
#
#    else:
#        blob_ell_l[idx] = size[size>0].mean()*1e-2
#        # Determine mean radial blob velocity in SOL
#        blob_vel_l[idx] = trail.get_velocity_com(rz_array)[good_pos_idx, 0].mean()

plt.title('$n_e / n_G = %3.2f,\, %d\, \\mathrm{events}$' % (n, blob_count) , fontsize=30)
#plt.title( '$n_e / n_G = %3.2f,\, %d\, \\mathrm{events}$' % ( n, blob_count ) )

corr_vr_ell = ( ( blob_ell - blob_ell.mean() ) * ( blob_vel - blob_vel.mean() ) ).mean() / \
    ( blob_vel.std() * blob_ell.std() )


good_size = (blob_ell > 0 ) & ( blob_ell < 0.99*blob_ell.max())

blob_ell_star = blob_ell[good_size] / ell_star
blob_vel_star = blob_vel[good_size] / v_star
    

#plt.plot( blob_ell_star, blob_vel_star, 'bo', \
#    label = '$V/V_* = %3.2f \\pm %3.2f$' % \
#    ( blob_vel_star.mean(), blob_vel_star.std() ) )
plt.plot( blob_ell_star, blob_vel_star, 'bo', \
    label = '$\\mathbf{V/V_* = %3.2f \\pm %3.2f}$' % \
    ( blob_vel_star.mean(), blob_vel_star.std() ) )

plt.ylim((0.0, 1.0))
leg = plt.legend(loc = 'upper right')

frame = leg.get_frame() 
frame.set_alpha(0.4)

plt.show()