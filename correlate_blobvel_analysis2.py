#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt

shotlist = [ 1120217010, 1120217014, 1120217015, 1120217020, 1120217021 ]

df_010 = np.load('1120217010/1120217010_blobvels_correlation.npz')
df_014 = np.load('1120217014/1120217014_blobvels_correlation.npz')
df_015 = np.load('1120217015/1120217015_blobvels_correlation.npz')
df_020 = np.load('1120217020/1120217020_blobvels_correlation.npz')
df_021 = np.load('1120217021/1120217021_blobvels_correlation.npz')

blobvel_010 = df_010['blob_vel']
blobvel_014 = df_014['blob_vel']
blobvel_015 = df_015['blob_vel']
blobvel_020 = df_020['blob_vel']
blobvel_021 = df_021['blob_vel']

blobvel_list = [blobvel_010, blobvel_014, blobvel_015, blobvel_020, blobvel_021]

for blobvel, shot in zip( blobvel_list, shotlist ):
    nonzero_idx = blobvel!=0
    mean_first = blobvel[nonzero_idx].mean()
    close_vels = np.abs( blobvel - mean_first) < blobvel.std()
    mean_second = blobvel[ close_vels ].mean()
    
    plt.figure()
    plt.plot( blobvel[close_vels].T, '.' )
    plt.title('# %d, mean = %5.2f m/s' % ( shot, mean_second ) )
    plt.xlabel('Blob event no.')
    plt.ylabel('Radial velocity / m/s')
    
    print '# %d, mean radial velocity: %5.2f' % ( shot, mean_second )
    
plt.show()
