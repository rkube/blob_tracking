#!/opt/local/bin/python
# -*- Encoding: UTF-8 -*-

"""
========
geometry
========

Convert blobtrails from pixel quantities to geometrical quantities
"""

import numpy as np

def com(array, xx=None, yy=None):
    """
    Return the center of mass position of the array x_com and y_com
    x_com = int(x * array) / int(array),
    y_com = int(y * array) / int(array)
    If xx and yy are not specified, use x = 0,1,...,np.shape(array)[0],
    and y = 0, 1, ..., np.shape(array)[1]

    Returns
        x_com, y_com
    """
    array = array.astype('float')
    if (xx is None and yy is None):
        xx, yy = np.meshgrid(np.arange(0, array.shape[0], 1.0),
                             np.arange(0, array.shape[1], 1.0))
        # If dx and dy are not specified, assume a regularly spaced grid
        return ((xx * array).sum() / array.sum(),
                (yy * array).sum() / array.sum())
    else:
        # Compute the increments in x and y
        dx = np.zeros_like(xx)
        dy = np.zeros_like(yy)
        dx[:, :-1] = xx[:, 1:] - xx[:, :-1]
        dx[:, -1] = dx[:, -2]

        dy[:-1, :] = yy[1:, :] - yy[:-1, :]
        dy[-2, :] = dy[-1, :]
        # Surface element
        dA = np.abs(dx) * np.abs(dy)
        return ((xx * array * dA).sum() / (array * dA).sum(),
                (yy * array * dA).sum() / (array * dA))


def com_rz(array, RR, zz):
    """
    Return the center of mass position on the irregulary spaced RR, zz array:
    R_com = int ( R * n * dA ) / int ( n * dA ), along second dimension
    z_com = int ( z * n * dA ) / int ( n * dA ), along first dimension
    """
    array = array.astype("float")
    dR, dz = np.zeros_like(RR), np.zeros_like(zz)
    dR[:, :-1] = RR[:, 1:] - RR[:, :-1]
    dR[:, -1] = dR[:, -2]
    dz[:-1, :] = zz[1:, :] - zz[:-1, :]
    dz[-1, :] = dz[:, -2]

    dA = np.abs(dR) * np.abs(dz)

    # COM along second dimension, COM along first dimension
    return((RR * array * dA).sum() / (array * dA).sum(),
           (zz * array * dA).sum() / (array * dz).sum())
    # return np.sum( RR * array * dR * dz ) / np.sum (array * dR * dz ) ,\
    #     np.sum( zz * array * dR * dz ) / np.sum (array * dR * dz )

#def find_sol_pixels(shotnr, frame_info=None, rz_array=None,
#                    datadir='/Users/ralph/source/blob_tracking/test_data'):
def find_sol_pixels(s):
    """
    Returns the indices of the pixels in between the separatrix and the LCFS.

    s:    Processed separatrix data from IDL
          i.e. s = readsav('%s/separatrix.sav' % (datadir), verbose=False)
          see /usr/local/cmod/codes/efit/idl/efit_rz2rmid.pro
              /home/terry/gpi/phantom/retrieve_phantom_RZ_array.pro
              /home/rkube/IDL/separatrix.pro,
    """


    gap_idx_mask = ((s['rmid'].reshape(64, 64) > s['rmid_sepx']) &
                    (s['rmid'].reshape(64, 64) < s['rmid_lim']))

    return np.argwhere(gap_idx_mask)


def find_sol_mask(shotnr, frame_info=None, rz_array=None,
                  datadir='/Users/ralph/source/blob_tracking/test_data'):
    """
    Returns a mask for the pixels in between the separatrix and the LCFS.
    """
    s = readsav('%s/separatrix.sav' % (datadir), verbose=False)

    return ((s['rmid'].reshape(64, 64) > s['rmid_sepx']) &
            (s['rmid'].reshape(64, 64) < s['rmid_lim']))


def blob_in_sol(trail, good_domain, logger=None):
    """
    Returns a bool array of the indices, in which the COM of a blob is
    in the SOL (good_domain)
    """

    try:
        # Consider only the positions, where the blob is in the good domain
        blob_pos = trail.get_trail_com()
        good_pos_idx = np.array([i in good_domain for i in
                                 blob_pos.round().astype('int').tolist()])

    except:
        good_pos_idx = np.ones_like(trail.get_tau())
        if (logger is not None):
            logger.info('This should not happen. Ignoring trigger domain')

    good_pos_idx = good_pos_idx[:-1]
    return good_pos_idx


def trail_com(trail, rz_array):
    """
    Return the position of the blob COM. Either in pixel or in (R,Z) coordinates if rz_array
    is passed.
    """

    if ( rz_array == None ):
        return self.xycom

    return rz_array[self.xycom[:,0].astype('int'), self.xycom[:,1].astype('int'), :]


def get_trail_max(self, rz_array = None):
    """
    Return the position of the blob maximum. Either in pixel or in (R,Z) coordinates if rz_array
    is passed.
    """
    if ( rz_array == None ):
        return self.xymax

    # Remember xycom[:,1] is the radial (X) index which corresponds to R
    return rz_array[self.xymax[:,0].astype('int'), self.xymax[:,1].astype('int'), :]


def velocity_max(trail, rz_array=None, dt=2.5e-6):
    """
    Return the velocity of the blob maximum. 
    Either in pixel / frame of m/s when rz_array is given
    
    Input:
        trail:    blobtrail
        rz_array: ndarray, axis0=radial, axis1=poloidal, axis2=(R,Z) in cm

    Output:
        vmax:     ndarray, either in m/s or px/frame
    
    """
    assert (np.size(trail.tau) > 1), 'Cannot compute blob velocity with only one frame recognized'
    xymax = trail.get_xymax().astype('int')

    try:
        vmax = 1e-2 * (rz_array[xymax[1:, 0], xymax[1:, 1], :]  - rz_array[xymax[:-1, 0], xymax[:-1, 1], :]) / dt
    except TypeError:
        vmax = xymax[1:, :] - xymax[:-1, :]

    return (vmax)

def velocity_com(trail, rz_array=None, dt=2.5e-6):
    """
    Return the velocity of the blob COM. Either in pixel / frame of m/s when rz_array is given

    Input:
        trail:    blobtrail
        rz_array: ndarray, axis0=radial, axis1=poloidal, axis2=(R,Z) in cm

    Output:
        vmax:     ndarray, either in m/s or px/frame
    """

    assert (np.size(trail.tau) > 1), 'Cannot compute blob velocity with only one frame recognized'
    xycom = trail.get_xycom().astype('int')

    try:
        vmax = 1e-2 * (rz_array[xycom[1:, 0], xycom[1:, 1], :]  - rz_array[xycom[:-1, 0], xycom[:-1, 1], :]) / dt
    except TypeError:
        vmax = xycom[1:, :] - xycom[:-1, :]

    return (vmax)


# End of file geometry.py
