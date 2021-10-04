# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 11:55:37 2021
Scattering function
@author: CHTUNG
"""
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
tStart = time.time()
chain01.chain()
tEnd = time.time()

#%% Calculate scattering function
DP = chain01.DP
chain_box = chain01.box
box_size = np.max(chain_box[1,:]-chain_box[0,:],axis=0)
n_grid = 512
grid_size = (box_size+1)/n_grid
Cc_relative = chain01.Cc.T-chain_box[0,:] # relative position of WL-chain in the box
bead_coord = np.floor(Cc_relative/grid_size).astype('int')

rho_r = np.zeros((n_grid,n_grid,n_grid))
for i in range(DP):
    rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1

rho_lmn = np.fft.fftn(rho_r)
S_lmn = np.absolute(rho_lmn)**2/DP

grid = np.meshgrid(np.arange(n_grid),np.arange(n_grid),np.arange(n_grid))
q_grid = np.sqrt(grid[0]**2+grid[1]**2+grid[2]**2)
index_q = np.floor(q_grid)
q_max = 2*np.pi*n_grid/box_size
S_q = np.zeros(int(n_grid/2))
qq = np.arange(int(n_grid/2))+1

for iq in range(int(n_grid/2)):
    vq = 4*np.pi*qq[iq]
    S_q[iq] = np.sum(S_lmn[index_q==iq])/vq

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6, 6),dpi=192)
ax = fig.add_subplot()
ax.plot(qq,S_q)
ax.set_xscale('log')
ax.set_yscale('log')
plt.show()
