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

    n_q = 50
    # qq = np.zeros(n_q)
    # qq = (np.logspace(-2,-0,n_q))
    qq = (np.linspace(1/n_q,1,n_q)/10)
    S_q = np.zeros(n_q)
    S_q_2D = np.zeros((n_q*2+1,n_q*2+1,3))
    
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
        #chain01.scatter_direct_block(qq,n_merge=1)
        chain01.scatter_direct_aniso(qq,n_merge=2)
        S_q = S_q + chain01.S_q
        S_q_2D = S_q_2D + chain01.S_q_2D
        tEnd = time.time()
        #print("\'scatter\' cost %f sec" % (tEnd - tStart))
    
    tEnd_loop = time.time()
    #print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))
    
    qq = chain01.qq
    qq_2D = chain01.qq_2D   
    S_q = S_q/n_chain
    S_q_2D = S_q_2D/n_chain
    
    # chain01.close()
    chain01.plot(axeslabel='on')
    
    filename = 'scatter_chain_2D_{:03f}_{:03f}_RB.mat'.format(qb,list_shear[j])
    mdic = {'S_q_2D':S_q_2D, 'S_q':S_q, 'qq':qq}
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
    


# #%% Calculate scattering function
# N = chain01.N
# chain_box = chain01.box
# approx_1D = 0

# #box_size = np.max(chain_box[1,:]-chain_box[0,:],axis=0)
# box_size = N
# n_grid = 256
# grid_size = (box_size)/n_grid
# Cc_relative = chain01.Cc.T-chain_box[0,:] # relative position of WL-chain in the box
# bead_coord = np.floor(Cc_relative/grid_size).astype('int')

# if approx_1D==1:
#     # density in real space
#     rho_rx = np.zeros(n_grid)
#     rho_ry = np.zeros(n_grid)
#     rho_rz = np.zeros(n_grid)
    
#     for i in range(N):
#         rho_rx[bead_coord[i,0]] += 1
#         rho_ry[bead_coord[i,0]] += 1
#         rho_rz[bead_coord[i,0]] += 1
    
#     # FFT and calculate scattering function
#     rho_qx = np.fft.fft(rho_rx)
#     rho_qy = np.fft.fft(rho_ry)
#     rho_qz = np.fft.fft(rho_rz)
#     S_q_x = np.absolute(rho_qx)**2/N
#     S_q_y = np.absolute(rho_qy)**2/N
#     S_q_z = np.absolute(rho_qz)**2/N
#     S_q_ave = (S_q_x + S_q_y + S_q_z)/3
    
#     # radial average
#     grid_coord = np.meshgrid(np.arange(n_grid))
#     dq_grid = 2*np.pi/(box_size)
#     q_grid = np.sqrt(grid_coord[0]**2)*dq_grid
    
#     dq = dq_grid
#     nq = int(np.floor(dq_grid/dq*n_grid/2))
#     qq = (np.arange(nq)+0.5)*dq    
#     index_q = np.floor(q_grid/dq) # q to the origin
    
#     S_q = np.zeros(int(nq))
    
#     for iq in range(int(nq)):
#         S_q[iq] = np.sum(S_q_ave[index_q==iq])/N
#         S_q[iq] = np.average(S_q_ave[index_q==iq])/N
    
# else:
#     # density in real space
#     rho_r = np.zeros((n_grid,n_grid,n_grid))
    
#     for i in range(N):
#         rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1
    
#     # FFT and calculate scattering function
#     rho_q = np.fft.fftn(rho_r)
#     S_q_lmn = np.absolute(rho_q)**2/N
    
#     # radial average
#     grid_coord = np.meshgrid(np.arange(n_grid),np.arange(n_grid),np.arange(n_grid))
#     dq_grid = 2*np.pi/(box_size)
#     q_grid = np.sqrt(grid_coord[0]**2+grid_coord[1]**2+grid_coord[2]**2)*dq_grid
    
#     dq = dq_grid
#     nq = int(np.floor(dq_grid/dq*n_grid/2))
#     qq = (np.arange(nq)+0.5)*dq
#     index_q = np.floor(q_grid/dq) # q to the origin
    
#     S_q = np.zeros(int(nq))
    
#     for iq in range(int(nq)):
#         #vq = 4*np.pi*(iq+0.5)**2/8
#         #S_q[iq] = np.sum(S_q_lmn[index_q==iq])/vq/N
#         S_q[iq] = np.average(S_q_lmn[index_q==iq])/N


