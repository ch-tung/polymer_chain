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
filename = 'scatter_chain_block_0.mat'
scatter_dict = loadmat(filename)
qq_max = 1
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']


n_mat = 10
S_q = np.zeros((np.shape(S_q_0)[0],np.shape(S_q_0)[1]*n_mat))
p = np.zeros((np.shape(p_0)[0],np.shape(p_0)[1]*n_mat))
for i in range(n_mat):
    filename = 'scatter_chain_block_{:d}.mat'.format(i)
    scatter_dict = loadmat(filename)
    S_q[:,i*np.shape(S_q_0)[1]:(i+1)*np.shape(S_q_0)[1]] = scatter_dict['S_q'][qq<qq_max,:]
    p[:,i*np.shape(p_0)[1]:(i+1)*np.shape(p_0)[1]] = scatter_dict['p']
    
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

qq = qq[qq<qq_max]

#%% PCA
index_p_ra = (p[0] == set_ra[0])
index_p_a2 = (p[1] == set_a2[0])
index_p_f = (p[2] == set_f[0])
# index_p = index_p_a2
index_p = np.arange(len(p[1])) # all datapoints

index_edge_ra_a2_0 = (p[0] == set_ra[0])&(p[1] == set_a2[0])
index_edge_ra_a2_1 = (p[0] == set_ra[0])&(p[1] == set_a2[len(set_a2)-1])
index_edge_ra_a2_2 = (p[0] == set_ra[len(set_ra)-1])&(p[1] == set_a2[0])
index_edge_ra_a2_3 = (p[0] == set_ra[len(set_ra)-1])&(p[1] == set_a2[len(set_a2)-1])
index_edge_ra_f_0 = (p[0] == set_ra[0])&(p[2] == set_f[0])
index_edge_ra_f_1 = (p[0] == set_ra[0])&(p[2] == set_f[len(set_f)-1])
index_edge_ra_f_2 = (p[0] == set_ra[len(set_ra)-1])&(p[2] == set_f[0])
index_edge_ra_f_3 = (p[0] == set_ra[len(set_ra)-1])&(p[2] == set_f[len(set_f)-1])
index_edge_f_a2_0 = (p[2] == set_f[0])&(p[1] == set_a2[0])
index_edge_f_a2_1 = (p[2] == set_f[0])&(p[1] == set_a2[len(set_a2)-1])
index_edge_f_a2_2 = (p[2] == set_f[len(set_f)-1])&(p[1] == set_a2[0])
index_edge_f_a2_3 = (p[2] == set_f[len(set_f)-1])&(p[1] == set_a2[len(set_a2)-1])

# property to be presented 0
pc_0 = p[0][index_p]
pc_0 = (pc_0-min(pc_0))/(max(pc_0)-min(pc_0))

# property to be presented 1
pc_1 = p[1][index_p]
pc_1 = (pc_1-min(pc_1))/(max(pc_1)-min(pc_1))

# property to be presented 2
pc_2 = p[2][index_p]
pc_2 = (pc_2-min(pc_2))/(max(pc_2)-min(pc_2))

c = plt.get_cmap('viridis')(pc_0)
c[:,0] = pc_0*1
c[:,1] = pc_1*0
c[:,2] = pc_2*0
c[:,3] = np.ones(c[:,0].shape)

# SVD
F = S_q
F = (F[:,:].T*qq).T
F = np.log(F)
# F = F - np.mean(F,axis=0)
U, S, Vh = np.linalg.svd(F)

score_F = np.matmul(F.T,U)

#%% plot
plt.close('all')

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(projection='3d')
# ax.plot(score_F[index_edge_0,0], score_F[index_edge_0,1], score_F[index_edge_0,2],'-k',lw=4)
# ax.plot(score_F[index_edge_1,0], score_F[index_edge_1,1], score_F[index_edge_1,2],'-k',lw=4)
# ax.plot(score_F[index_edge_2,0], score_F[index_edge_2,1], score_F[index_edge_2,2],'-k',lw=4)
# ax.plot(score_F[index_edge_3,0], score_F[index_edge_3,1], score_F[index_edge_3,2],'-k',lw=4)

# connect points with the same ra
for i in range(len(set_ra)):
    index_edge_ra = (p[0] == set_ra[i])&index_p
    ax.plot(score_F[index_edge_ra,0], score_F[index_edge_ra,1], score_F[index_edge_ra,2],'-r')
    
for i in range(len(set_f)):
    index_edge_f = (p[2] == set_f[i])&index_p
    ax.plot(score_F[index_edge_f,0], score_F[index_edge_f,1], score_F[index_edge_f,2],'-b')

ax.scatter(score_F[index_p,0], score_F[index_p,1], score_F[index_p,2], 
            'o',
            s=10,
            alpha=1,
            lw=2,
            facecolors=c,
            edgecolors=c)
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