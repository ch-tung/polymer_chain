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
import threading

#%% test
# backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
N_backbone = 5000

# a1 = N_backbone/np.array([1.25,10,80,640])
a2 = N_backbone/np.array([1.25,10,80,640])
ra = 2**np.arange(8)
f = 1/2**(np.arange(8))
p = np.meshgrid(ra,a2,f)
p_ra = p[0].flatten()
p_a2 = p[1].flatten()
p_f = p[2].flatten()
n_p = len(p_f)

def scattering_loop(n_q,n_chain,chain01):
    """
    generate polymer chains of the same parameters and calculate scattering function
    Parameters
    ----------
    n_q : int
        number of q points.
    n_chain : int
        number of chains used to calculate the averaged S(Q).
    chain01 : class
        WLM class

    Returns
    -------
    qq : TYPE
        q.
    S_q_j : TYPE
        S(Q).

    """
    qq = 2*np.pi/(np.logspace(1,5,n_q))
    S_q_j = np.zeros(n_q)
    for i in np.arange(n_chain):

        # tStart = time.time()
        chain01.apply_SA = 0
        chain01.d_exc = chain01.a*0.1*2
        chain01.chain_block()
        
        # tStart = time.time()
        # chain01.scatter_grid(n_grid=n_q*2)
        # chain01.scatter_grid_direct(n_q=len(qq),n_grid=256,box_size=np.max(chain_box[1,:]-chain_box[0,:])+1) 
        chain01.scatter_direct(qq,n_merge=2)
        S_q_j = S_q_j + chain01.S_q
        # tEnd = time.time()
        # print("\'scatter\' cost %f sec" % (tEnd - tStart))
        
    S_q_j = S_q_j/n_chain
        
    return qq, S_q_j

def job(j):
    # parameters
    # Chain stiffness
    a_backbone = np.array([p_a2[j]*p_ra[j],p_a2[j]])
    
    # block ratio
    f = p_f[j]
	
    # Unit persistence
    lambda_backbone = 1
    
    # call class
    chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
    chain01.f = f
    
    n_q = 64 
    n_chain = 100

    # tStart_loop = time.time()
    qq, S_q_j = scattering_loop(n_q,n_chain,chain01)
    # tEnd_loop = time.time()
    # print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))
    
    S_q[:,j] = S_q_j
    
S_q = np.zeros((64,n_p))
n_j = n_p

class MyThread(threading.Thread):
    def __init__(self, j):
        threading.Thread.__init__(self)
        self.j = j
    
    def run(self):
        job(self.j)

threads = []
tStart = time.time()

for j in range(n_j):
    
    # job(j)
    # threads.append(threading.Thread(target = job, args = (j,)))
    threads.append(MyThread(j))
    threads[j].start()
    
for j in range(n_j):
    threads[j].join()

qq = 2*np.pi/(np.logspace(1,5,64))

tEnd = time.time()
print("it cost %f sec" % (tEnd - tStart))
from scipy.io import savemat
filename = 'scatter_chain_block.mat'
mdic = {'S_q':S_q, 'p':np.array([p_ra,p_a2,p_f]), 'qq':qq}
savemat(filename, mdic)