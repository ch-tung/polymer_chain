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
DP_backbone = 10000

# Chain stiffness
a_backbone = 1e2

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(DP_backbone,a_backbone,lambda_backbone,unit_C)
chain01.d_exc = 0.5

n_q = 64
qq = np.zeros(64)
S_q = np.zeros(64)

n_chain = 10
for i in range(n_chain):
    
    tStart = time.time()
    #chain01.apply_SA = 0
    chain01.chain()
    tEnd = time.time()
    print("\'chain\' cost %f sec" % (tEnd - tStart))
    
    tStart = time.time()
    chain01.scatter(n_grid=n_q*2)
    S_q = S_q + chain01.S_q
    tEnd = time.time()
    print("\'scatter\' cost %f sec" % (tEnd - tStart))

qq = chain01.qq    
S_q = S_q/n_chain


import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6, 6),dpi=192)
ax = fig.add_subplot()
ax.plot(qq,S_q)
ax.set_xscale('log')
ax.set_yscale('log')
plt.show()

# #%% Calculate scattering function
# DP = chain01.DP
# chain_box = chain01.box
# approx_1D = 0

# #box_size = np.max(chain_box[1,:]-chain_box[0,:],axis=0)
# box_size = DP
# n_grid = 256
# grid_size = (box_size)/n_grid
# Cc_relative = chain01.Cc.T-chain_box[0,:] # relative position of WL-chain in the box
# bead_coord = np.floor(Cc_relative/grid_size).astype('int')

# if approx_1D==1:
#     # density in real space
#     rho_rx = np.zeros(n_grid)
#     rho_ry = np.zeros(n_grid)
#     rho_rz = np.zeros(n_grid)
    
#     for i in range(DP):
#         rho_rx[bead_coord[i,0]] += 1
#         rho_ry[bead_coord[i,0]] += 1
#         rho_rz[bead_coord[i,0]] += 1
    
#     # FFT and calculate scattering function
#     rho_qx = np.fft.fft(rho_rx)
#     rho_qy = np.fft.fft(rho_ry)
#     rho_qz = np.fft.fft(rho_rz)
#     S_q_x = np.absolute(rho_qx)**2/DP
#     S_q_y = np.absolute(rho_qy)**2/DP
#     S_q_z = np.absolute(rho_qz)**2/DP
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
#         S_q[iq] = np.sum(S_q_ave[index_q==iq])/DP
#         S_q[iq] = np.average(S_q_ave[index_q==iq])/DP
    
# else:
#     # density in real space
#     rho_r = np.zeros((n_grid,n_grid,n_grid))
    
#     for i in range(DP):
#         rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1
    
#     # FFT and calculate scattering function
#     rho_q = np.fft.fftn(rho_r)
#     S_q_lmn = np.absolute(rho_q)**2/DP
    
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
#         #S_q[iq] = np.sum(S_q_lmn[index_q==iq])/vq/DP
#         S_q[iq] = np.average(S_q_lmn[index_q==iq])/DP


