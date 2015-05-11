#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

"""
==============
test1
==============

.. codeauthor:: Ralph Kube <ralphkube@gmail.com>

Test 2d gauss fitting on exemplary data
"""


import numpy as np
import matplotlib.pyplot as plt
import gaussfitter as gf
from scipy.ndimage import gaussian_filter
from functions import gauss2d, com
import time

datafile = np.load('1100803015/frames_mean.npz')
frames_mean = datafile['frames_mean']

tic = time.clock()
datafile = np.load('1100803015/frames.npz')
toc = time.clock()

latexout = True
smoothing   = True

print 'Loading frames. Time elapsed: ', toc - tic
# Read frames and flip horizontally so that large x is radially outwards
#frames_data = datafile['frames_data'][:,:,::-1]

x = np.arange(0, 64,)
y = np.arange(0, 64)
xx, yy = np.meshgrid(x, y)

frame0 = 30332
nframes = 12
frame_range = np.arange(frame0, frame0 + nframes)

A0_t = np.zeros(nframes)
x_t = np.zeros(nframes)
y_t = np.zeros(nframes)
sx_t = np.zeros(nframes)
sy_t = np.zeros(nframes)
xc_t = np.zeros(nframes)
yc_t = np.zeros(nframes)

i = 0

if ( latexout == True ):
    print '\\begin{tabular}{c|c|c|c|c|c|c}'
    print 'Frameno. & $A_0$ & $x_0$ & $y_0$ & $\\sigma_x$   & $\\sigma_y$   & $\\theta$ \\\\'
for frameno in frame_range:

    # To subtract the mean, we have to convert from uint16 to int32 (signed!).
    # do this in-place
    test_frame_u = datafile['frames_data'][frameno, :, ::-1] 
    test_frame = test_frame_u.view('int16')
    test_frame[:,:] = test_frame_u
    test_frame = test_frame - frames_mean

    if ( smoothing == True ):
        test_frame = gaussian_filter(test_frame, (1,1) )

    plt.figure()
    plt.title('t = %d' % frameno)
    plt.contourf(xx, yy, test_frame, 16)
    plt.colorbar()

    A0 = np.max(test_frame)
    j, ix = np.unravel_index( np.argmax(test_frame), np.shape(test_frame) )
#    print 'Max: %5.2f at (%d,%d)' % (A0, jy, ix)

    p0 = (0., A0, ix, jy, 1., 1., 0. )
    tic = time.clock()
    result = gf.gaussfit( test_frame, params = p0, return_all = 1 )
    toc = time.clock()
#    print 'Fitting. Time elapsed: ', toc - tic

    A0_t[i] = result[0][1]
    x_t[i]= result[0][2]
    y_t[i] = result[0][3]
    sx_t[i] = result[0][4]
    sy_t[i] = result[0][5]

    g2 = gauss2d( xx, yy, result[0] )
    xc_t[i], yc_t[i] = com( xx, yy, g2 )

    plt.contour( xx, yy, g2 )
    if ( latexout == True ):
        print '$%d$ & $%4.2f \\pm %4.2f$    & $%4.2f \\pm %4.2f$    & $%4.2f \\pm %4.2f$    & $%3.2f \\pm %3.2f$    & $%3.2f \\pm %3.2f$    & $%4.1f \\pm %4.1f$    \\\\'\
            % (frameno, result[0][1], result[1][1], result[0][2], result[1][2], \
                result[0][3], result[1][3], result[0][4], result[1][4], \
                result[0][5], result[1][5], result[0][6], result[1][6] )


    i = i + 1 
    test_frame, test_frame_u = 0.0, 0.0

print '\\end{tabular}'
 


plt.figure()
plt.title('Amplitude')
plt.plot(frame_range, A0_t)

plt.figure()
plt.title('Position')
plt.plot(frame_range, x_t, label='x')
plt.plot(frame_range, y_t, label='y')
plt.plot(frame_range, xc_t, label='XC')
plt.plot(frame_range, yc_t, label='YC')

plt.legend()

plt.figure()
plt.title('Width')
plt.plot(frame_range, sx_t, label='$\\sigma_x$')
plt.plot(frame_range, sy_t, label='$\\sigma_y$')
plt.legend()


plt.show()

