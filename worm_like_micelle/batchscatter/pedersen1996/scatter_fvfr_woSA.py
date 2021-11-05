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

n_a = 4
a = N_backbone/np.array([1.25,10,80,640])

def scattering_loop(n_q,n_chain,chain01):
    qq = 2*np.pi/(np.logspace(1,5,n_q))
    S_q_j = np.zeros(n_q)
    for i in np.arange(n_chain):

        # tStart = time.time()
        chain01.apply_SA = 0
        chain01.d_exc = chain01.a*0.1*2
        # chain01.chain()
        chain01.chain_fix_val_free_rot()
        # chain01.ring(n_harmonics=40,sigma=10)
        # chain01.ring_q()
        # tEnd = time.time()
        # print("\'chain\' cost %f sec" % (tEnd - tStart))
        
        # N = chain01.N
        # chain_box = chain01.box
        
        # tStart = time.time()
        # chain01.scatter_grid(n_grid=n_q*2)
        # chain01.scatter_grid_direct(n_q=len(qq),n_grid=256,box_size=np.max(chain_box[1,:]-chain_box[0,:])+1) 
        chain01.scatter_direct(qq,n_merge=1)
        S_q_j = S_q_j + chain01.S_q
        # tEnd = time.time()
        # print("\'scatter\' cost %f sec" % (tEnd - tStart))
        
    S_q_j = S_q_j/n_chain
        
    return qq, S_q_j

def job(j):
    # Chain stiffness
    a_backbone = a[j]
    # a_backbone = 2e3
	
    # Unit persistence
    lambda_backbone = 1
    
    # call class
    chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
    #chain01.d_exc = 1
    
    n_q = 64 
    n_chain = 1000

    # tStart_loop = time.time()
    qq, S_q_j = scattering_loop(n_q,n_chain,chain01)
    # tEnd_loop = time.time()
    # print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))
    
    S_q[:,j] = S_q_j

S_q = np.zeros((64,n_a))
n_j = n_a

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
filename = 'scatter_chain_fvfr_woSA.mat'
mdic = {'S_q':S_q, 'a':a, 'qq':qq}
savemat(filename, mdic)