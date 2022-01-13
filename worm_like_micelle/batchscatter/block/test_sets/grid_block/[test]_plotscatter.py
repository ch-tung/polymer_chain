# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 22:15:38 2021
PCA on S(Q)
@author: CHTUNG
"""

import numpy as np
import numpy.matlib
from scipy.io import loadmat
# from scipy.io import savemat
# import matplotlib
import matplotlib.pyplot as plt
# from scipy import interpolate

#%% load
# load mat files
# grep shape
filename = 'scatter_chain_block_1_g.mat'
scatter_dict = loadmat(filename)
qq_max = 2
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']

filename = '../grid_block_SA/scatter_chain_block_1_g.mat'
scatter_dict = loadmat(filename)
qq_max = 2
qq = scatter_dict['qq'][0,:]
S_q_0_SA  = scatter_dict['S_q'][qq<qq_max,:]
p_0_SA = scatter_dict['p']

# fix the deviation high q limit
import Sk
q = qq
L = 1000
b = L*2
S_q_rod = Sk.S_rod(q,L)

S_q_rod_bead = loadmat('scatter_rod.mat')['S_q']

delta_S_q = S_q_rod/S_q_rod_bead

# n_homo = 100
# S_q_homo = np.zeros((len(q),n_homo))
# for i in range(n_homo):
#     S_q_homo[:,i] = Sk.Sk(q,L,L/(0.1+(i)*100/n_homo))

d_th = 100
S_q_original = S_q_0*1
# S_q = (S_q.T*delta_S_q).T
# for i in range(S_q.shape[1]):
#     chi = np.exp(-(q*d_th)**(-5))
#     S_q_i = S_q[:,i]
#     # S_q[:,i] = chi*S_q_rod + (1-chi)*S_q_i
#     S_q_i[S_q_i<S_q_rod] = (chi*S_q_rod)[S_q_i<S_q_rod] + ((1-chi)*S_q_i)[S_q_i<S_q_rod]
#     # S_q_i[S_q_i<S_q_rod] = S_q_rod[S_q_i<S_q_rod]
#     S_q[:,i] = S_q_i

# S_q_homo_original = S_q_homo*1
# S_q_homo = (S_q_homo.T*delta_S_q).T
# for i in range(S_q_homo.shape[1]):
#     chi = np.exp(-(q*d_th)**(-5))
#     S_q_i = S_q_homo[:,i]
#     # S_q_homo[:,i] = chi*S_q_rod + (1-chi)*S_q_homo_i
#     S_q_i[S_q_i<S_q_rod] = (chi*S_q_rod)[S_q_i<S_q_rod] + ((1-chi)*S_q_i)[S_q_i<S_q_rod]
#     S_q_homo[:,i] = S_q_i

qq = qq[qq<qq_max]

#%% plot S_q
fig = plt.figure(figsize=(6, 6))
ax_SQ = fig.add_subplot()

ax_SQ.plot(qq.T*L,S_q_rod,'--k')
for i in range(16):
    ax_SQ.plot((10**(2*i)*np.array([1e-4, 1e1]))**-(1/4)*L,np.array([1e-4, 1e1]),
        '--',color='#C0C0C0',linewidth=0.5)
    ax_SQ.plot((10**(i)*np.array([1e-4, 1e1]))**-(1/2)*L,np.array([1e-4, 1e1]),
            ':',color='#C0C0C0',linewidth=0.5)
    # ax_SQ.plot((10**(i)*np.array([1e-4, 1e1]))**-(0.588)*L,np.array([1e-4, 1e1]),
    #         ':',color='#C0C0C0',linewidth=0.5)
    ax_SQ.plot((10**(i)*np.array([1e-4, 1e1]))**-(1)*L,np.array([1e-4, 1e1]),
            '-.',color='#C0C0C0',linewidth=0.5)

for i in range(1):
    ax_SQ.plot(qq.T*L,S_q_0[:,90],'-')
for i in range(1):
    ax_SQ.plot(qq.T*L,S_q_0_SA[:,90],'--')

ax_SQ.set_xscale('log')
ax_SQ.set_yscale('log')
ax_SQ.set_xlabel('QL')
ax_SQ.set_ylabel('S(Q)')
ax_SQ.set_xlim([1e-4*L, 1e0*L])
ax_SQ.set_ylim([0.5e-2, 2e0])

