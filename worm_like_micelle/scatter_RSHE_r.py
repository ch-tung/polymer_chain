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
N_backbone = 1000

# Chain stiffness
# a_backbone = np.array([1e2,1e1])
a_backbone = 1e1

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
chain01.d_exc = 1

n_q = 65
# qq = np.zeros(n_q)
qq = (np.logspace(-4,-0,n_q))
n_r = 31
# qq = np.zeros(n_q)
rr = (np.logspace(0,3,n_r))
# qq = (np.linspace(1/n_q,1,n_q)/10)
S_q = np.zeros(n_q)
S_q_lm = np.zeros((int(n_q),6))

g_r = np.zeros(n_r)
g_r_lm = np.zeros((int(n_r),6))

n_chain = 100
tStart_loop = time.time()
for i in range(n_chain):
    print('{:03d}/{:03d}'.format(i+1,n_chain))
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
    qb = 0.005
    shear = 0.00
    chain01.kappa = -np.log(qb)
    chain01.epsilon = shear
    chain01.grid = 'SC'
    chain01.chain_grid_shear()
    tEnd = time.time()
    print("\'chain\' cost %f sec" % (tEnd - tStart))
    
    N = chain01.N
    chain_box = chain01.box
    
    tStart = time.time()
    # chain01.scatter_grid(n_grid=n_q*2)
    # chain01.scatter_grid_direct(n_q=len(qq),n_grid=256,box_size=np.max(chain_box[1,:]-chain_box[0,:])+1) 
    # chain01.scatter_direct_pw(qq,n_merge=1)
    # chain01.scatter_direct_block(qq,n_merge=1)
    # chain01.scatter_direct_aniso(qq,n_merge=2)
    chain01.scatter_direct_RSHE(qq,rr,n_merge=1,calculate_g_r=1)
    g_r = g_r + chain01.g_r
    g_r_lm = g_r_lm + chain01.g_r_lm
    tEnd = time.time()
    print("\'scatter\' cost %f sec" % (tEnd - tStart))

tEnd_loop = time.time()
print("\'loop\' cost %f sec" % (tEnd_loop - tStart_loop))

rr = chain01.rr
g_r = g_r/n_chain
g_r_lm = g_r_lm/n_chain

# chain01.close()
chain01.plot(axeslabel='on')

#%%
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6, 6),dpi=96)
ax = fig.add_subplot()

from scipy.io import loadmat

ax.plot(rr,g_r)
ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim([np.sqrt(np.min(qq)*np.max(qq))*10**-2,np.sqrt(np.min(qq)*np.max(qq))*10**2])
# ax.set_ylim([1e-3, 2e0])
ax.grid(True,which='major')

#%%
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8, 6),dpi=96)
fig.tight_layout()
plt.rcParams['font.size'] = '16'
ax = fig.add_subplot()
gs00 = ax.plot(rr,g_r_lm[:,0],'-',color='#000000',label='$(0,0)$')
s2n2 = ax.plot(rr,g_r_lm[:,1],'-',color='#FF0000',label='$(2,-2)$')
s2n1 = ax.plot(rr,g_r_lm[:,2],'-',color='#00FF00',label='$(2,-1)$')
s20 = ax.plot(rr,g_r_lm[:,3],'-',color='#0000FF',label='$(2,0)$')
s2p1 = ax.plot(rr,g_r_lm[:,4],'-',color='#00FFFF',label='$(2,1)$')
s2p2 = ax.plot(rr,g_r_lm[:,5],'-',color='#FF00FF',label='$(2,2)$')

ax.legend()

ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim([np.sqrt(np.min(qq)*np.max(qq))*10**-2,np.sqrt(np.min(qq)*np.max(qq))*10**2])
# ax.set_ylim([-2e-1, 1.1e0])
ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$g_l^{m}(r)$')
ax.grid(True,which='major')

# from scipy.io import savemat
# filename = 'scatter_chain_aniso_{:03f}_{:03f}.mat'.format(qb,shear)
# mdic = {'S_q_lm':S_q_lm, 'S_q':S_q, 'qq':qq}
# savemat(filename, mdic)
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

