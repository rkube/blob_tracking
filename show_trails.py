#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import logging
import cPickle
#import matplotlib as mpl
#mpl.use('AGG')
import matplotlib.pyplot as plt
from load_mdsframes import load_mdsframes
import helper_functions
import blob_tracking
import blob_statistics

frame0 = 0
save = False
shotnr = 1091216036

# Store information about density scan with loggin facility
logger = 0.
logger = logging.getLogger('density_scan_091216')
logger.setLevel(logging.DEBUG)
file_handle = logging.FileHandler('logs/ds091216_trails.log')
file_handle.setLevel(logging.DEBUG)
formatter = logging.Formatter( '%(message)s' )
#file_handle.setFormatter(formatter)
logger.addHandler(file_handle)

logger.info('Start logging!!!')

print 'Processing shot# %d' % shotnr
logger.info('Processing shot %d' % shotnr)    

datafile = np.load('%d/%d_testframes.npz' % (shotnr, shotnr) )
frame_info = datafile['frame_info']

picklefile = open('%d/%d_trails.pkl' % (shotnr, shotnr), 'rb' )
trails = cPickle.load(picklefile)
logger.info('Loaded blob trails from file')
picklefile.close()


good_domain = helper_functions.find_sol_pixels(shotnr, frame_info)
#print type(good_domain)
good_domain = np.array(good_domain)[:,:].tolist()
#print good_domain[:,:].tolist()
print good_domain
good_domain_arr = np.array(good_domain)

num_events = len(trails)

for idx, trail in enumerate(trails):
    try:
        # Consider only the positions, where the blob is in the good domain
        blob_pos = trail.get_trail_com()
        good_pos_idx = np.array([i in good_domain for i in blob_pos.round().astype('int').tolist()])   
        
#        print [i in good_domain for i in blob_pos.round().astype('int').tolist()]
        
    except:
        good_pos_idx = np.ones_like(trail.get_tau())
        if ( logger != None ):
            logger.info('This should not happen. Ignoring trigger domain for blob %d' % idx)

    good_pos_idx = good_pos_idx[:-1]
#    print good_pos_idx, good_pos_idx.sum()
    # Skip this blob if it is not detected between LCFS and LS
    if ( good_pos_idx.sum() < 4 ):
        if ( logger != None ):
            logger.info('Blob %d/%d was not detected in the SOL, rejecting' % (idx, num_events) )
        else:
            print 'Blob %d/%d was not detected in the SOL, rejecting' % (idx, num_events)
        continue

#    print (trail.get_trail_com())[:,0], (trail.get_trail_com())[:,1]


    plt.figure()
    plt.title('%d' % good_pos_idx.sum() )
    plt.plot(good_domain_arr[:,1], good_domain_arr[:,0], 'ko')
    plt.plot( trail.get_trail_com()[:,1], trail.get_trail_com()[:,0] ,'x')


    plt.show()
