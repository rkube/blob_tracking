#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
from load_mdsframes import load_mdsframes
import matplotlib.pyplot as plt
from scipy.stats import kurtosis, skew


"""
Compute and plot
1. <I>_RMS (R,Z)    RMS as a function of R,Z
2. <I>_t (R,Z)      Time avg. as a function of R,Z
3. S(I) (R,Z)       Skewness as a function of R,Z
4. K(I) (R,Z)       Kurtosis as a function of R,Z
"""
def comp_sk(shotnr, nframes):

    try:
        print 'Loading data from %d/%d_skew_kurt.npz' % ( shotnr, shotnr )
        datafile = np.load('%d/%d_skew_kurt.npz' %(shotnr, shotnr) )
        shotnr = datafile['shotnr']
        i_skew = datafile['i_skew']
        i_kurt = datafile['i_kurt']
        frame0 = datafile['frame0']
    except IOError, KeyError :
        print 'Could not find data in specified file. Computing from phantom frames'
        frame0 = 0 
    
    #    datafile = np.load('%d/%d_frames.npz' % (shotnr, shotnr) )
    #    if ( datafile['transposed'] == False ):
    #        frames = datafile['frames_data'][frame0:frame0+nframes, :, ::-1]
    #    else:
    #        frames = datafile['frames_data'][frame0 : frame0+nframes, :, :]
    
        frame_range = np.arange(frame0, frame0+nframes)
    
        print 'Computing skewness'
        i_skew = skew(frames, axis=0)
        print 'Computing kurtosis'
        i_kurt = kurtosis(frames, axis=0, fisher=True)
    
    plt.figure()
    plt.title('#%d, Skewness' % shotnr)
    plt.xlabel('R / px')
    plt.ylabel('Z / px')
    plt.contourf(i_skew)
    plt.colorbar()
    
    plt.figure()
    plt.title('#%d, Skewness z=32' % shotnr)
    plt.xlabel('R / px')
    plt.plot(i_skew[32,:])
    
    
    plt.figure()
    plt.title('#%d, Kurtosis' % shotnr)
    plt.contourf(i_kurt)
    plt.xlabel('R / px')
    plt.ylabel('Z / px')
    plt.colorbar()
    
    plt.figure()
    plt.title('#%d, Kurtosis z=32' % shotnr)
    plt.xlabel('R / px')
    plt.plot(i_kurt[32,:])
    
    plt.show()

#np.savez('%d/%d_skew_kurt.npz' % (shotnr, shotnr) , i_skew=i_skew, i_kurt=i_kurt, shotnr=shotnr, frame0=frame0, transposed=True)


def plot_sk(shotnr, skew, kurt):

    plt.figure()
    plt.title('#%d, Skewness' % shotnr)
    plt.xlabel('R / px')
    plt.ylabel('Z / px')
    plt.contourf(i_skew)
    plt.colorbar()
    
    plt.figure()
    plt.title('#%d, Skewness z=32' % shotnr)
    plt.xlabel('R / px')
    plt.plot(i_skew[32,:])
    
    
    plt.figure()
    plt.title('#%d, Kurtosis' % shotnr)
    plt.contourf(i_kurt)
    plt.xlabel('R / px')
    plt.ylabel('Z / px')
    plt.colorbar()
    
    plt.figure()
    plt.title('#%d, Kurtosis z=32' % shotnr)
    plt.xlabel('R / px')
    plt.plot(i_kurt[32,:])
    
    plt.show()


