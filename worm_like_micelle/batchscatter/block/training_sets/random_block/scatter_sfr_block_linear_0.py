# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 11:55:37 2021
Scattering function
@author: CHTUNG
"""
# import sys
# sys.modules[__name__].__dict__.clear()

import numpy as np
import numpy.matlib
import time
from WLM import WLChain
import threading
# import multiprocessing

#%% test
# backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
N_backbone = 1000

# load parameters
from scipy.io import loadmat
filename_stats = 'stats_block.mat'
filename_parameters = 'parameters_block.mat'
stats_dict = loadmat(filename_stats)
parameters_dict = loadmat(filename_parameters)

parameters = parameters_dict['parameters']
stats = stats_dict['statistics']

set_stats0 = sorted(set(stats[:,0]))
set_stats1 = sorted(set(stats[:,1]))
set_stats2 = sorted(set(stats[:,2]))

i_s0 = 0
n_set = 10
n_sample = 500
n_total = n_set*n_sample;
index_s0 = np.arange(i_s0*n_sample,(i_s0+1)*n_sample)
p_ra = parameters[index_s0,0]
p_a2 = (N_backbone*parameters[index_s0,1])
p_f = parameters[index_s0,2]
# p = np.meshgrid(ra,a2,f)
# p_ra = p[0].flatten()
# p_a2 = p[1].flatten()
# p_f = p[2].flatten()
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
    qq = (np.logspace(-4,0,65))
    S_q_j = np.zeros(n_q)
    for i in np.arange(n_chain):

        # tStart = time.time()
        chain01.apply_SA = 0
        chain01.d_exc = chain01.a*0.1*2
        chain01.chain_block()
        
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
    # print('This is Process: ', j)
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
    
    n_q = 65 
    n_chain = 1000

    # tStart_loop = time.time()
    qq, S_q_j = scattering_loop(n_q,n_chain,chain01)
    # tEnd_loop = time.time()
    # print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))
    
    S_q[:,j] = S_q_j
    
#%% run
    
S_q = np.zeros((65,n_p))
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

qq = (np.logspace(-4,0,65))

tEnd = time.time()
print("it cost %f sec" % (tEnd - tStart))
from scipy.io import savemat
filename = 'scatter_chain_block_{:d}_r.mat'.format(i_s0)
mdic = {'S_q':S_q, 'p':np.array([p_ra,p_a2,p_f]), 'qq':qq}
savemat(filename, mdic)