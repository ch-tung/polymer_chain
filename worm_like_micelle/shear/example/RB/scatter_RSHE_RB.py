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
from scipy.io import savemat

#%% test
# backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
N_backbone = 1000

# Chain stiffness
# a_backbone = np.array([1e2,1e1])
a_backbone = 1e1

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
chain01.d_exc = 1


qb = 0.0001
list_shear = [0.00,0.005,0.01,0.02,0.04]

def job(j):
    
    n_q = 65
    # qq = np.zeros(n_q)
    qq = (np.logspace(-4,-0,n_q))
    # qq = (np.linspace(1/n_q,1,n_q)/10)
    S_q = np.zeros(n_q)
    S_q_2D = np.zeros((n_q*2+1,n_q*2+1,3))
    S_q_lm = np.zeros((n_q,6))
    
    n_chain = 1000
    tStart_loop = time.time()
    for i in np.arange(n_chain):
        #print('{:03d}/{:03d}'.format(i+1,n_chain))
        tStart = time.time()
        chain01.apply_SA = 0
        chain01.d_exc = chain01.a*0.1*2
        chain01.f = 0.5
        # chain01.chain_block()
        # chain01.chain()
        # chain01.chain_fix_val_free_rot()
        #chain01.ring(n_harmonics=40,sigma=10)
        #chain01.ring_q()
        
        # chain_grid method
        chain01.d_exc = 1
        chain01.kappa = -np.log(qb)
        chain01.epsilon = list_shear[j]
        chain01.grid = 'RB'
        chain01.chain_grid_shear()
        tEnd = time.time()
        #print("\'chain\' cost %f sec" % (tEnd - tStart))
        
        N = chain01.N
        chain_box = chain01.box
        
        tStart = time.time()
        # chain01.scatter_grid(n_grid=n_q*2)
        # chain01.scatter_grid_direct(n_q=len(qq),n_grid=256,box_size=np.max(chain_box[1,:]-chain_box[0,:])+1) 
        # chain01.scatter_direct_pw(qq,n_merge=1)
        # chain01.scatter_direct_block(qq,n_merge=1)
        # chain01.scatter_direct_aniso(qq,n_merge=2)
        chain01.scatter_direct_RSHE(qq,n_merge=1)
        S_q = S_q + chain01.S_q
        S_q_lm = S_q_lm + chain01.S_q_lm
        # S_q_2D = S_q_2D + chain01.S_q_2D
        tEnd = time.time()
        #print("\'scatter\' cost %f sec" % (tEnd - tStart))
    
    tEnd_loop = time.time()
    #print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))
    
    qq = chain01.qq
    # qq_2D = chain01.qq_2D   
    S_q = S_q/n_chain
    S_q_lm = S_q_lm/n_chain
    S_q_2D = S_q_2D/n_chain
    
    filename = 'scatter_chain_aniso_{:03f}_{:03f}_RB.mat'.format(qb,list_shear[j])
    mdic = {'S_q_lm':S_q_lm, 'S_q':S_q, 'qq':qq}
    savemat(filename, mdic)

class MyThread(threading.Thread):
    def __init__(self, j):
        threading.Thread.__init__(self)
        self.j = j
    
    def run(self):
        job(self.j)

threads = []
tStart = time.time()

for j, shear in enumerate(list_shear):
    
    threads.append(MyThread(j))
    threads[j].start()
    
for j, shear in enumerate(list_shear):
    threads[j].join()

# #%%
# yy, xx = np.meshgrid(qq_2D,qq_2D)
# fig_2D, ax_2D = plt.subplots(1,3, figsize=[12, 4])
# fig_2D.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.4, hspace=None)
# labels = [[r'$Q_y$',r'$Q_z$'],[r'$Q_x$',r'$Q_z$'],[r'$Q_x$',r'$Q_y$']]

# for i in range(3):
#     ax_2D[i].pcolor(xx, yy, S_q_2D[:,:,i],shading='auto',vmin=0.0, vmax=1.0)
#     # ax_2D[i].pcolor(S_q_2D[:,:,i],shading='auto')
#     ax_2D[i].set_xlabel(labels[i][0])
#     ax_2D[i].set_ylabel(labels[i][1])
#     ax_2D[i].set_aspect('equal', 'box')
# plt.show()

