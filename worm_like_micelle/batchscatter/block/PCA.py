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
filename = 'scatter_chain_block_3.mat'
scatter_dict = loadmat(filename)
S_q  = scatter_dict['S_q']
p = scatter_dict['p']
qq = scatter_dict['qq']
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

#%% SVD
F = (S_q.T*qq).T
F = np.log(S_q)
# F = F - np.mean(F,axis=0)
F = F
U, S, Vh = np.linalg.svd(F)

score_F = np.matmul(F.T,U)

#%% plot
plt.close('all')
index_p_a2 = (p[1] == set_a2[0])
index_p_ra = (p[0] == set_ra[7])
index_p_f = (p[2] == set_f[6])
# index_p = index_p_ra
index_p = np.arange(len(p[1]))
pc = np.log(p[0][index_p])
# pc = p[2][index_p]

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(projection='3d')
ax.scatter(score_F[index_p,0], score_F[index_p,1], score_F[index_p,2], 
           c=pc,
           s=40,
           alpha=1,
           edgecolors = [0,0,0])
ax.view_init(elev=25, azim=-135)
ax.set_xlabel('SVD[0]')
ax.set_ylabel('SVD[1]')
ax.set_zlabel('SVD[2]')

plt.show()

#%% plot basis
fig = plt.figure(figsize=(6, 6))
ax_basis = fig.add_subplot()

n_basis = 4
x = np.linspace(0.0, 1.0, n_basis)
color = plt.get_cmap('viridis')(x)

for i in range(n_basis):
    ax_basis.plot(qq.T,U[:,i], color = color[i])

ax_basis.set_xscale('log')
ax_basis.set_xlabel('Q')
ax_basis.set_ylabel('score')

#%% plot variance
fig = plt.figure(figsize=(6, 6))
ax_var = fig.add_subplot()
ax_var.plot(np.arange(len(S)),S)

ax_var.set_yscale('log')
ax_var.set_xlabel('rank')
ax_var.set_ylabel(r'$\Sigma$')