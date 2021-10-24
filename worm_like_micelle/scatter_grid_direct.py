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

def scatter(Cc, n_grid=256, box_size=1e4):
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

    #box_size = np.max(chain_box[1,:]-chain_box[0,:],axis=0)
    #box_size = N
    grid_size = (box_size)/n_grid
    Cc_relative = Cc.T-chain_box[0,:] # relative position of WL-chain in the box
    bead_coord = np.floor(Cc_relative/grid_size).astype('int')

    # density in real space
    rho_r = np.zeros((n_grid,n_grid,n_grid))
    
    for i in range(N):
        rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1
    
    # FFT and calculate scattering function
    rho_q = np.fft.fftn(rho_r)
    S_q_lmn = np.absolute(rho_q)**2/N
    
    # radial average
    grid_coord = np.meshgrid(np.arange(n_grid),np.arange(n_grid),np.arange(n_grid))
    dq_grid = 2*np.pi/(box_size)
    q_grid = np.sqrt(grid_coord[0]**2+grid_coord[1]**2+grid_coord[2]**2)*dq_grid
    
    dq = dq_grid
    nq = int(np.floor(dq_grid/dq*n_grid/2))
    qq = (np.arange(nq)+0.5)*dq
    index_q = np.floor(q_grid/dq) # q to the origin
    
    S_q = np.zeros(int(nq))
    
    for iq in range(int(nq)):
        #vq = 4*np.pi*(iq+0.5)**2/8
        #S_q[iq] = np.sum(S_q_lmn[index_q==iq])/vq/N
        S_q[iq] = np.average(S_q_lmn[index_q==iq])/N
            
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
a_backbone = 2e2

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
chain01.d_exc = 1

n_q = 128
qq = np.zeros(n_q)
S_q = np.zeros(n_q)

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

S_q, qq = scatter(Cc)

tEnd = time.time()
print("\'scatter\' cost %f sec" % (tEnd - tStart))

tEnd_loop = time.time()
print("\'loop\' cost %f sec" % (tEnd - tStart))

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
