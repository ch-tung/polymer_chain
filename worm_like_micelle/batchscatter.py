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
S_q = np.zeros((64,10))
n_j = 10

tStart_loopj = time.time()
for j in range(n_j):
    # Chain stiffness
    a_backbone = n_j**((j+2)/2)
    
    # Unit persistence
    lambda_backbone = 1
    
    # call class
    chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
    chain01.d_exc = 1
    
    n_q = 64
    S_q_j = np.zeros(n_q)
    qq = 2*np.pi/(np.logspace(1,5,n_q))
    
    n_chain = 1
    tStart_loop = time.time()
    for i in range(n_chain):
    
        tStart = time.time()
        chain01.apply_SA = 1
        chain01.chain()
        #chain01.ring(n_harmonics=40,sigma=10)
        #chain01.ring_q()
        tEnd = time.time()
        print("\'chain\' cost %f sec" % (tEnd - tStart))
        
        N = chain01.N
        chain_box = chain01.box
        
        tStart = time.time()
        # chain01.scatter_grid(n_grid=n_q*2)
        # chain01.scatter_grid_direct(n_q=len(qq),n_grid=256,box_size=np.max(chain_box[1,:]-chain_box[0,:])+1) 
        chain01.scatter_direct(qq,n_merge=4)
        S_q_j = S_q_j + chain01.S_q
        tEnd = time.time()
        print("\'scatter\' cost %f sec" % (tEnd - tStart))
    
    tEnd_loop = time.time()
    print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))
    
    qq = chain01.qq    
    S_q_j = S_q_j/n_chain
    a[j] = a_backbone
    S_q[:,j] = S_q_j
	
tEnd_loopj = time.time()
print("it cost %f sec" % (tEnd_loopj - tStart_loopj))

#from scipy.io import savemat
#filename = 'scatter_chain_prstnc.mat'
#mdic = {'S_q':S_q, 'a':a, 'qq':qq}
#savemat(filename, mdic)