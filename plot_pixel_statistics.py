#!/usr/bin/env python
# -*- Encoding: UTF-8 -*-

# import MDSplus as mds
import numpy as np
import matplotlib.pyplot as plt

"""
Plot time series, histogram, and DFT of a given pixel
"""


def gauss2d(A, x, y, x0, y0, sigma_x, sigma_y, theta):
    a = 0.5 * ((np.cos(theta) / sigma_x) ** 2.0
               + (np.sin(theta) / sigma_y) ** 2.0)
    b = -0.25 * (np.sin(2.0 * theta) *
                 (sigma_x ** -2.0 + sigma_y ** -2.0))
    c = 0.5 * ((np.sin(theta) / sigma_x) ** 2. +
               (np.cos(theta) / sigma_y) ** 2)

    return A * np.exp(-a * (x - x0) ** 2 -
                      2.0 * b * (x - x0) * (y - y0)
                      - c * (y - y0) ** 2)


# myTree = mds.Tree('spectroscopy', 1100803015  )
# node2 = myTree.getNode('GPI.PHANTOM.FRAMES')
# frames = node2.getData()
# frames_data = frames.data()

shotnr = 1100803015
# datadir = '/Volumes/CMOD/phantom/'
datadir = '/Users/ralph/source/blob_tracking/'
df_fname = '%s/%10d/%10d_frames.npz' % (datadir, shotnr, shotnr)

with np.load(df_fname) as df:
    frames = df['frames_normalized_rms'][:, :, ::-1]
    df.close()

dt = 2.5e-6
tb = np.arange(0, frames.shape[0]) * dt
ix, jy = 52, 32

x = np.arange(0., 1., 1./64.)
y = np.arange(0., 1., 1./64.)
xx, yy = np.meshgrid(x, y)

A_l = np.zeros(50000)
x0_l = np.zeros(50000)
y0_l = np.zeros(50000)
sx_l = np.zeros(50000)
sy_l = np.zeros(50000)
t_l = np.zeros(50000)

pixel_str = '#%10d: jy=%02d, ix=%02d' % (shotnr, jy, ix)

# Compute Hisogram and Spectrum
t_ft = np.fft.rfft(frames[:, jy, ix])
freq = np.fft.fftfreq(frames.shape[0], dt)[:frames.shape[0] / 2]

hist, edges = np.histogram(frames[:, jy, ix])
mid = edges[:-1] + 0.5 * np.diff(edges).mean()

fig_ts = plt.figure()
ax_ts = fig_ts.add_subplot(111)
ax_ts.plot(tb, frames[:, jy, ix])
ax_ts.set_xlabel('time / microseconds')
ax_ts.set_ylabel('I / a.u.')
ax_ts.set_title(pixel_str)


fig_ft = plt.figure()
ax_ft = fig_ft.add_subplot(111)
ax_ft.plot(freq, t_ft * t_ft.conj())
ax_ft.set_xlabel('Frequency / Hz')
ax_ft.set_ylabel('|I(omega)|^2')
ax_ts.set_title(pixel_str)


fig_hs = plt.figure()
ax_hs = fig_hs.add_subplot(111)
ax_hs.plot(mid, hist, '.')
ax_hs.set_xlabel('I / a.u.')
ax_hs.set_ylabel('PDF(I)')
ax_hs.set_title(pixel_str)


A_l[23002] = 222
x0_l[23002] = 46. / 64.
y0_l[23002] = 1. / 64.
sx_l[23002] = 0.1
sy_l[23002] = 0.1
t_l[23002] = 0.0


for t in np.arange(23002, 23003):
    plt.figure()
    plt.title('t = %5d' % t)
    plt.contourf(xx, yy, frames[t, :, :])
    plt.contour(xx, yy, gauss2d(A_l[t], xx, yy, x0_l[t], x0_l[t],
                                sx_l[t], sy_l[t], t_l[t]), 8)
    # plt.plot(frames_data[t, jy, :])
    print frames[t, jy, ix]
    plt.colorbar()


plt.show()

# End of file analyze.py
