# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 11:55:37 2021
Scattering function
@author: CHTUNG
"""
import sys
sys.modules[__name__].__dict__.clear()

import numpy as np
import numpy.matlib
import time
from WLM import WLChain

#%% test
# backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
N_backbone = 10000

a = np.zeros(10)
S_q = np.zeros((128,10))
for i in range(10):
    # Chain stiffness
    a_backbone = 10**((i+1)/2)
    
    # Unit persistence
    lambda_backbone = 1
    
    # call class
    chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
    chain01.d_exc = 1
    
    n_q = 128
    qq = np.zeros(n_q)
    S_q = np.zeros(n_q)
    
    n_chain = 100
    tStart_loop = time.time()
    for i in range(n_chain):
    
        tStart = time.time()
        #chain01.apply_SA = 0
        #chain01.chain()
        #chain01.ring(n_harmonics=40,sigma=10)
        chain01.ring_q()
        tEnd = time.time()
        print("\'chain\' cost %f sec" % (tEnd - tStart))
        
        tStart = time.time()
        chain01.scatter(n_grid=n_q*2)
        S_q = S_q + chain01.S_q
        tEnd = time.time()
        print("\'scatter\' cost %f sec" % (tEnd - tStart))
    
    tEnd_loop = time.time()
    print("\'loop\' cost %f sec" % (tEnd - tStart))
    
    qq = chain01.qq    
    S_q_i = S_q/n_chain
    a[i] = a_backbone
    S_q[:,i] = S_q_i

from scipy.io import savemat
filename = 'scatter_chain.mat'
mdic = {'S_q':S_q, 'a':a, 'qq':qq}
savemat(filename, mdic)