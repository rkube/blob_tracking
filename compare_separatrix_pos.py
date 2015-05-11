#!/opt/local/bin/python
#-*- Encoding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import idlsave

"""
Plot the separatrix positions for the discharges on 100803 and 120217
"""



shotlist = [ 1120711020, 1120711022, 1120711024, 1120711027 ]
    
ne_ng = [ 0.2, 0.31, 0.32, 0.40 ]
    
    
seplist = np.zeros( len(shotlist) )
    
for shot, idx in zip(shotlist, np.arange( len(shotlist) ) ):
    sep = idlsave.read('%d/%d_separatrix.sav' % (shot, shot) , verbose = False)
    seplist[idx] = sep.rmid_sepx * 100
    
plt.plot(ne_ng, seplist, 'bs', label='120711')
#plt.xlabel('Shot')
#plt.ylabel('Separatrix at Zmid / cm')

#plt.figure()
#plt.xlabel('Shot')
plt.xlabel('$n_e / n_G$')
plt.ylabel('Separatrix at Zmid / cm')

plt.legend()
plt.show()