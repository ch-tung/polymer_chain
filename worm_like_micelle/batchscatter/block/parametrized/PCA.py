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

# load parameters
from scipy.io import loadmat
filename_stats = 'stats_block.mat'
filename_parameters = 'parameters_block.mat'
stats_dict = loadmat(filename_stats)
parameters_dict = loadmat(filename_parameters)

parameters = parameters_dict['parameters']
stats = stats_dict['statistics']

set_stats0 = sorted(set(stats[:,0]))
set_stats1 = sorted(set(stats[:,1]))
set_stats2 = sorted(set(stats[:,2]))

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
index_p_ra = (p[0] != set_ra[0])
index_p_a2 = (p[1] == set_a2[0])
index_p_f = (p[2] == set_f[0])
index_p = index_p_ra
index_p = np.arange(len(p[1])) # all datapoints

# index_edge_ra_a2_0 = (p[0] == set_ra[0])&(p[1] == set_a2[0])
# index_edge_ra_a2_1 = (p[0] == set_ra[0])&(p[1] == set_a2[len(set_a2)-1])
# index_edge_ra_a2_2 = (p[0] == set_ra[len(set_ra)-1])&(p[1] == set_a2[0])
# index_edge_ra_a2_3 = (p[0] == set_ra[len(set_ra)-1])&(p[1] == set_a2[len(set_a2)-1])
# index_edge_ra_f_0 = (p[0] == set_ra[0])&(p[2] == set_f[0])
# index_edge_ra_f_1 = (p[0] == set_ra[0])&(p[2] == set_f[len(set_f)-1])
# index_edge_ra_f_2 = (p[0] == set_ra[len(set_ra)-1])&(p[2] == set_f[0])
# index_edge_ra_f_3 = (p[0] == set_ra[len(set_ra)-1])&(p[2] == set_f[len(set_f)-1])
# index_edge_f_a2_0 = (p[2] == set_f[0])&(p[1] == set_a2[0])
# index_edge_f_a2_1 = (p[2] == set_f[0])&(p[1] == set_a2[len(set_a2)-1])
# index_edge_f_a2_2 = (p[2] == set_f[len(set_f)-1])&(p[1] == set_a2[0])
# index_edge_f_a2_3 = (p[2] == set_f[len(set_f)-1])&(p[1] == set_a2[len(set_a2)-1])

# property to be presented 0
p_0 = p[0][index_p]
pc_0 = (p_0-np.nanmin(p_0))/(np.nanmax(p_0)-np.nanmin(p_0))

# property to be presented 1
p_1 = p[1][index_p]/1000
pc_1 = (p_1-np.nanmin(p_1))/(np.nanmax(p_1)-np.nanmin(p_1))

# property to be presented 2
p_2 = p[2][index_p]
pc_2 = (p_2-np.nanmin(p_2))/(np.nanmax(p_2)-np.nanmin(p_2))

# first moment
pm_1 = np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2


# second moment
pm_2 = (np.log(p_1)**2*(1-p_2) + np.log(p_0*p_1)**2*p_2 - 
    (np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2)**2)


# third moment
pm_3 = (p_2*(np.log(p_0*p_1)-pm_1)**3 + (1-p_2)*(np.log(p_1)-pm_1)**3)/np.sqrt(pm_2)**3


pm_c_1 = (pm_1-np.nanmin(pm_1))/(np.nanmax(pm_1)-np.nanmin(pm_1))
pm_c_2 = (pm_2-np.nanmin(pm_2))/(np.nanmax(pm_2)-np.nanmin(pm_2))
pm_c_3 = (pm_3-np.nanmin(pm_3[np.isfinite(pm_3)]))/(np.nanmax(pm_3[np.isfinite(pm_3)])-np.nanmin(pm_3[np.isfinite(pm_3)]))

# SVD
F = S_q
F = (F[:,:].T*qq).T
F = np.log(F)
F = F - np.mean(F,axis=0)
# F = np.gradient(F,axis=0)
# F = np.gradient(F,axis=0)
U, S, Vh = np.linalg.svd(F)

score_F = np.matmul(F.T,U)

#%% plot
plt.close('all')

#%% plot SVD
c = plt.get_cmap('viridis')(pm_c_3)
# c[:,0] = pc_0*0
# c[:,1] = pc_1*0
# c[:,2] = pc_2*1
# c[:,3] = np.ones(c[:,0].shape)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(projection='3d')
# ax.plot(score_F[index_edge_0,0], score_F[index_edge_0,1], score_F[index_edge_0,2],'-k',lw=4)
# ax.plot(score_F[index_edge_1,0], score_F[index_edge_1,1], score_F[index_edge_1,2],'-k',lw=4)
# ax.plot(score_F[index_edge_2,0], score_F[index_edge_2,1], score_F[index_edge_2,2],'-k',lw=4)
# ax.plot(score_F[index_edge_3,0], score_F[index_edge_3,1], score_F[index_edge_3,2],'-k',lw=4)

# connect points with the same ra
# for i in range(len(set_ra)):
#     index_edge_ra = (p[0] == set_ra[i])&index_p
#     ax.plot(score_F[index_edge_ra,0], score_F[index_edge_ra,1], score_F[index_edge_ra,2],'-r')
    
# for i in range(len(set_f)):
#     index_edge_f = (p[2] == set_f[i])&index_p
#     ax.plot(score_F[index_edge_f,0], score_F[index_edge_f,1], score_F[index_edge_f,2],'-b')

ax.scatter(score_F[index_p,0], score_F[index_p,1], score_F[index_p,2], 
            'o',
            s=10,
            alpha=1,
            lw=2,
            facecolors=c,
            edgecolors=c)
ax.view_init(elev=25, azim=25)
ax.set_xlabel('SVD[0]')
ax.set_ylabel('SVD[1]')
ax.set_zlabel('SVD[2]')

plt.show()

#%% plot basis
fig = plt.figure(figsize=(6, 6))
ax_basis = fig.add_subplot()

n_basis = 3
x = np.linspace(0.0, 1.0, n_basis)
color = plt.get_cmap('viridis')(x)

for i in range(n_basis):
    ax_basis.plot(qq.T,U[:,i], color = color[i])

ax_basis.set_xscale('log')
ax_basis.set_xlabel('Q')
ax_basis.set_ylabel('score')

#%% plot score
score_c = score_F[index_p,1]
score_n = (score_c-np.nanmin(score_c))/(np.nanmax(score_c)-np.nanmin(score_c))
c_score = plt.get_cmap('viridis')(score_n)

#%% plot score in parameter space
# plt.close('all')

fig = plt.figure(figsize=(6, 6))
ax_s = fig.add_subplot(projection='3d')

# connect points with the same ra
ax_s.scatter(np.log(p_1), np.log(p_0), p_2, 
            'o',
            s=100,
            c = c_score,
            alpha=1,
            lw=0,
            edgecolors=c_score)
ax_s.view_init(elev=22.5, azim=-150)
ax_s.set_xlabel(r'$lnb_1$')
ax_s.set_ylabel(r'$lnr_b$')
ax_s.set_zlabel(r'$f$')

plt.show()

#%% plot score in moment space
# plt.close('all')

fig = plt.figure(figsize=(6, 6))
ax_s = fig.add_subplot(projection='3d')

# connect points with the same ra
ax_s.scatter(stats[index_p,0], stats[index_p,1], stats[index_p,2], 
            'o',
            s=100,
            c = c_score,
            alpha=1,
            lw=0,
            edgecolors=c_score)
ax_s.view_init(elev=20, azim=-150)
ax_s.set_xlabel(r'$\mu$')
ax_s.set_ylabel(r'$\sigma^{2}$')
ax_s.set_zlabel(r'$\gamma$')

plt.show()

#%% plot scree
fig = plt.figure(figsize=(6, 6))
ax_var = fig.add_subplot()
ax_var.plot(np.arange(len(S)),S)

ax_var.set_yscale('log')
ax_var.set_xlabel('rank')
ax_var.set_ylabel(r'$\Sigma$')
