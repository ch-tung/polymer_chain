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

def scatter_direct(Cc, n_grid=512, box_size=1e4):
    """
    Calculate scattering function.
    
    Args:
        n_grid: int
            number of grid points
        approx_1D: boolean
            1-D FFT for isotropic systems
    """
    
    # N = self.N
    # chain_box = self.box

    grid_size = (box_size)/n_grid
    Cc_relative = Cc.T-chain_box[0,:] # relative position of WL-chain in the box
    bead_coord = np.floor(Cc_relative/grid_size).astype('int')
    
    # density in real space
    rho_r = np.zeros((n_grid,n_grid,n_grid))
    
    for i in range(N):
        rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1
    
    list_rho_r_all = rho_r.reshape(-1)
    index_rho_r = np.nonzero(list_rho_r_all)[0]
    list_rho_r = list_rho_r_all[index_rho_r]
    coord_rho_r_x = np.floor_divide(index_rho_r, n_grid**2)
    coord_rho_r_y = np.floor_divide(index_rho_r-coord_rho_r_x*n_grid**2, n_grid)
    coord_rho_r_z = index_rho_r-coord_rho_r_x*n_grid**2-coord_rho_r_y*n_grid
    coord_rho_r = np.vstack((coord_rho_r_x,coord_rho_r_y,coord_rho_r_z))*box_size/n_grid
    
    # two-point correlation
    n_list = len(list_rho_r)
    r_jk = coord_rho_r.reshape(n_list,1,3) - coord_rho_r.reshape(1,n_list,3)
    d_jk = np.sqrt(np.sum(r_jk**2,axis=2))
    n_jk = np.outer(list_rho_r, list_rho_r)
    rho_jk = n_jk/np.sum(n_jk)
    
    # radial average
    # dq_grid = 2*np.pi/(box_size)
    # dq = dq_grid
    # nq = int(np.floor(dq_grid/dq*n_grid/2))
    # qq0 = np.arange(nq)+0.5
    # qq = qq0*dq
    qq0 = np.logspace(0,3,32)*2*np.pi/1e5
    nq = len(qq0)
    qq = qq0 
    
    S_q = np.zeros(int(nq))
    d_jk_list = d_jk[d_jk!=0]
    rho_jk_list = rho_jk[d_jk!=0]
    
    for iq in range(int(nq)):
        sinqr_qr = rho_jk_list*np.sin(qq0[iq]*d_jk_list)/(qq0[iq]*d_jk_list)
        S_q[iq] = np.sum(sinqr_qr[np.isnan(sinqr_qr)==0])
            
    #self.qq = qq
    #self.S_q = S_q
    return S_q, qq

#%% test
# backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
N_backbone = 10000

# Chain stiffness
a_backbone = 1e2

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
chain01.d_exc = 1
chain01.close()

# n_q = 128
# qq = np.zeros(n_q)
# S_q = np.zeros(n_q)

tStart_loop = time.time()

tStart = time.time()
#chain01.apply_SA = 0
chain01.chain()
tEnd = time.time()
print("\'chain\' cost %f sec" % (tEnd - tStart))

tStart = time.time()

N = chain01.N
chain_box = chain01.box
Cc = chain01.Cc

S_q, qq = scatter_direct(Cc,n_grid=512)

# n_grid=512
# box_size=1e4

# grid_size = (box_size)/n_grid
# Cc_relative = Cc.T-chain_box[0,:] # relative position of WL-chain in the box
# bead_coord = np.floor(Cc_relative/grid_size).astype('int')

# # density in real space
# rho_r = np.zeros((n_grid,n_grid,n_grid))

# for i in range(N):
#     rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1

# list_rho_r_all = rho_r.reshape(-1)
# index_rho_r = np.nonzero(list_rho_r_all)[0]
# list_rho_r = list_rho_r_all[index_rho_r]
# coord_rho_r_x = np.floor_divide(index_rho_r, n_grid**2)
# coord_rho_r_y = np.floor_divide(index_rho_r-coord_rho_r_x*n_grid**2, n_grid)
# coord_rho_r_z = index_rho_r-coord_rho_r_x*n_grid**2-coord_rho_r_y*n_grid
# coord_rho_r = np.vstack((coord_rho_r_x,coord_rho_r_y,coord_rho_r_z))*box_size/n_grid

# # two-point correlation
# n_list = len(list_rho_r)
# r_jk = coord_rho_r.reshape(n_list,1,3) - coord_rho_r.reshape(1,n_list,3)
# d_jk = np.sqrt(np.sum(r_jk**2,axis=2))
# n_jk = np.outer(list_rho_r, list_rho_r)
# rho_jk = n_jk/np.sum(n_jk)

# # radial average
# # dq_grid = 2*np.pi/(box_size)
# # dq = dq_grid
# # nq = int(np.floor(dq_grid/dq*n_grid/2))
# # qq0 = np.arange(nq)+0.5
# # qq = qq0*dq
# qq0 = np.logspace(0,3,32)*2*np.pi/1e5
# nq = len(qq0)
# qq = qq0 

# S_q = np.zeros(int(nq))
# d_jk_list = d_jk[d_jk!=0]
# rho_jk_list = rho_jk[d_jk!=0]

# for iq in range(int(nq)):
#     sinqr_qr = rho_jk_list*np.sin(qq0[iq]*d_jk_list)/(qq0[iq]*d_jk_list)
#     S_q[iq] = np.sum(sinqr_qr[np.isnan(sinqr_qr)==0])

tEnd = time.time()
print("\'scatter\' cost %f sec" % (tEnd - tStart))

tEnd_loop = time.time()
print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))

chain01.plot()

#%%
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6, 6),dpi=192)
ax = fig.add_subplot()

for i in range(5):
    ax.plot((10**(i+3)*np.array([1e-3, 1e1]))**-(1/2),np.array([1e-3, 1e1]),
            '--',color='#C0C0C0',linewidth=0.5)
    ax.plot((10**(i+3)*np.array([1e-3, 1e1]))**-(3/5),np.array([1e-3, 1e1]),
            ':',color='#C0C0C0',linewidth=0.5)

ax.plot(qq,S_q)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([np.sqrt(np.min(qq)*np.max(qq))*10**-1.5,np.sqrt(np.min(qq)*np.max(qq))*10**1.5])
ax.set_ylim([1e-3, 2e0])
ax.grid(True,which='major')
plt.show()
