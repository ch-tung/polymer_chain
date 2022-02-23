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
filename = 'scatter_chain_block_0_r.mat'
scatter_dict = loadmat(filename)
qq_max = 2
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']

# load parameters
from scipy.io import loadmat
# filename_stats = 'stats_block.mat'
filename_parameters = 'parameters_block.mat'
# stats_dict = loadmat(filename_stats)
parameters_dict = loadmat(filename_parameters)

parameters = parameters_dict['parameters']

# load stats
# filename_stats = 'stats_block.mat'
filename_parameters = 'parameters_block.mat'
# stats_dict = loadmat(filename_stats)
parameters_dict = loadmat(filename_parameters)

# stats = stats_dict['statistics']
# p_m1 = stats[:,0]
# p_m2 = stats[:,1]
# p_m3 = stats[:,2]

# set_m1 = sorted(set(p_m1))
# set_m2 = sorted(set(p_m2))
# set_m3 = sorted(set(p_m3))

n_mat = 10
S_q = np.zeros((np.shape(S_q_0)[0],np.shape(S_q_0)[1]*n_mat))
p = np.zeros((np.shape(p_0)[0],np.shape(p_0)[1]*n_mat))
n_set = 10
n_sample = 500
for i in range(n_mat):
    # index_s0 = stats[:,0]==set_m1[i]
    index_s0 = np.arange(i*n_sample,(i+1)*n_sample)
    filename = 'scatter_chain_block_{:d}_r.mat'.format(i)
    scatter_dict = loadmat(filename)
    S_q[:,index_s0] = scatter_dict['S_q'][qq<qq_max,:]
    p[:,index_s0] = scatter_dict['p']
    

# fix the deviation high q limit
import Sk
q = qq
L = 1000
b = L*2
S_q_rod = Sk.S_rod(q,L)

S_q_rod_bead = loadmat('scatter_rod.mat')['S_q']

delta_S_q = S_q_rod/S_q_rod_bead

n_homo = 100
S_q_homo = np.zeros((len(q),n_homo))
for i in range(n_homo):
    S_q_homo[:,i] = Sk.Sk(q,L,L/(0.1+(i)*100/n_homo))

d_th = 100
S_q_original = S_q*1
S_q = (S_q.T*delta_S_q).T
# for i in range(S_q.shape[1]):
#     chi = np.exp(-(q*d_th)**(-5))
#     S_q_i = S_q[:,i]
#     # S_q[:,i] = chi*S_q_rod + (1-chi)*S_q_i
#     S_q_i[S_q_i<S_q_rod] = (chi*S_q_rod)[S_q_i<S_q_rod] + ((1-chi)*S_q_i)[S_q_i<S_q_rod]
#     # S_q_i[S_q_i<S_q_rod] = S_q_rod[S_q_i<S_q_rod]
#     S_q[:,i] = S_q_i

S_q_homo_original = S_q_homo*1
# S_q_homo = (S_q_homo.T*delta_S_q).T
# for i in range(S_q_homo.shape[1]):
#     chi = np.exp(-(q*d_th)**(-5))
#     S_q_i = S_q_homo[:,i]
#     # S_q_homo[:,i] = chi*S_q_rod + (1-chi)*S_q_homo_i
#     S_q_i[S_q_i<S_q_rod] = (chi*S_q_rod)[S_q_i<S_q_rod] + ((1-chi)*S_q_i)[S_q_i<S_q_rod]
#     S_q_homo[:,i] = S_q_i
    
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

qq = qq[qq<qq_max]

#%% PCA
# property to be presented 0 (rb)
p_0 = p[0]
pc_0 = (p_0-np.nanmin(p_0))/(np.nanmax(p_0)-np.nanmin(p_0))

# property to be presented 1 (b1)
p_1 = p[1]/1000
pc_1 = (p_1-np.nanmin(p_1))/(np.nanmax(p_1)-np.nanmin(p_1))

# property to be presented 2 (f)
p_2 = p[2]
pc_2 = (p_2-np.nanmin(p_2))/(np.nanmax(p_2)-np.nanmin(p_2))

# first moment
pm_1 = np.log10(p_1)*(1-p_2) + np.log10(p_0*p_1)*p_2


# second moment
pm_2 = (np.log10(p_1)**2*(1-p_2) + np.log10(p_0*p_1)**2*p_2 - 
    (np.log10(p_1)*(1-p_2) + np.log10(p_0*p_1)*p_2)**2)


# third moment
pm_3 = (p_2*(np.log10(p_0*p_1)-pm_1)**3 + (1-p_2)*(np.log10(p_1)-pm_1)**3)/np.sqrt(pm_2)**3


pm_c_1 = (pm_1-np.nanmin(pm_1))/(np.nanmax(pm_1)-np.nanmin(pm_1))
pm_c_2 = (pm_2-np.nanmin(pm_2))/(np.nanmax(pm_2)-np.nanmin(pm_2))
pm_c_3 = (pm_3-np.nanmin(pm_3[np.isfinite(pm_3)]))/(np.nanmax(pm_3[np.isfinite(pm_3)])-np.nanmin(pm_3[np.isfinite(pm_3)]))

# index_p_m1 = (p_m1 == set_m1[10])
# index_p_m2 = (p_m2 == set_m2[10])
# index_p_m3 = (p_m3 == set_m3[10])
index_p_ra = (p[0] == set_ra[4])
index_p_a2 = (p[1] == set_a2[0])
index_p_f = (p[2] == set_f[4])
index_p_all = np.arange(len(p[1])) # all datapoints, bool
index_p = index_p_all[p[0][index_p_all]<100]

# SVD
def feature(S_q_in):
    F = S_q_in
    F = (F.T/S_q_rod.T).T
    # F = (F[:,:].T*qq).T
    F = np.log(F)
    # F = F - np.mean(F,axis=0)
    # F = np.gradient(F,axis=0)
    return F

F = feature(S_q)
F_rod = feature(S_q_rod)
F_homo = feature(S_q_homo)
U, S, Vh = np.linalg.svd(F[:,index_p])

score_F = np.matmul(F.T,U)
score_F_rod = np.matmul(F_rod.T,U)
score_F_homo = np.matmul(F_homo.T,U)

#%% plot
plt.close('all')

#%% plot SVD
# c[:,0] = pc_0*0
# c[:,1] = pc_1*0
# c[:,2] = pc_2*1
# c[:,3] = np.ones(c[:,0].shape)
c_plot_SVD = np.log(p_1[index_p])
x_plot_SVD = ((c_plot_SVD-np.nanmin(c_plot_SVD[np.isfinite(c_plot_SVD)])) /
          (np.nanmax(c_plot_SVD[np.isfinite(c_plot_SVD)])-np.nanmin(c_plot_SVD[np.isfinite(c_plot_SVD)])))

color_SVD = plt.get_cmap('viridis')(x_plot_SVD)

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
            facecolors=color_SVD,
            edgecolors=color_SVD)
ax.scatter(0, 0, 0, 
            'o',
            s=20,
            alpha=1,
            lw=2,
            facecolors='k',
            edgecolors='k')
ax.plot(score_F_homo[:,0], score_F_homo[:,1], score_F_homo[:,2], 
            '-r',
            lw=4)
ax.view_init(elev=36, azim=45)
ax.set_xlabel('SVD[0]')
ax.set_ylabel('SVD[1]')
ax.set_zlabel('SVD[2]')
d_score = np.max(score_F[index_p,0:3],axis = 0)-np.min(score_F[index_p,0:3],axis = 0)
ax.set_box_aspect([d_score[0],d_score[1],d_score[2]])

plt.show()

#%% plot basis
fig = plt.figure(figsize=(6, 6))
ax_basis = fig.add_subplot()

n_basis = 3
x = np.linspace(0.0, 1.0, n_basis)
color_basis = plt.get_cmap('viridis')(x)
for i in range(n_basis):
    ax_basis.plot(qq.T*L,U[:,i], color = color_basis[i])

ax_basis.set_xscale('log')
ax_basis.set_xlabel('QL')
ax_basis.set_ylabel('score')

#%% plot F

# index to plot
index_plot = index_p_a2&index_p_ra

# property to be presented in color
c_plot = (p_2[index_plot])
x_plot = ((c_plot-np.nanmin(c_plot[np.isfinite(c_plot)])) /
          (np.nanmax(c_plot[np.isfinite(c_plot)])-np.nanmin(c_plot[np.isfinite(c_plot)])))

color = plt.get_cmap('viridis')(x_plot)

fig = plt.figure(figsize=(6, 6))
ax_F = fig.add_subplot()

ax_F.plot(qq.T*L,0*qq.T,'--k')
for i in range(np.sum(index_plot)):
    ax_F.plot(qq.T*L,F[:,index_p_all[index_plot][i]], color = color[i])

ax_F.set_xscale('log')
ax_F.set_xlabel('QL')
ax_F.set_ylabel('F')
ax_F.set_xlim([1e-4*L, 1e0*L])

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

for i in range(np.sum(index_plot)):
    ax_SQ.plot(qq.T*L,S_q[:,index_p_all[index_plot][i]], color = color[i])

ax_SQ.set_xscale('log')
ax_SQ.set_yscale('log')
ax_SQ.set_xlabel('QL')
ax_SQ.set_ylabel('S(Q)')
ax_SQ.set_xlim([1e-4*L, 1e0*L])
ax_SQ.set_ylim([0.5e-2, 2e0])

#%% plot score
# score_c = score_F[index_p_all,1]
# score_n = (score_c-np.nanmin(score_c))/(np.nanmax(score_c)-np.nanmin(score_c))
# c_score = plt.get_cmap('viridis')(score_n)

# # plot score in parameter space
# # plt.close('all')

# fig = plt.figure(figsize=(6, 6))
# ax_s = fig.add_subplot(projection='3d')

# # connect points with the same ra
# ax_s.scatter(np.log(p_1), np.log(p_0), p_2, 
#             'o',
#             s=100,
#             c = c_score,
#             alpha=1,
#             lw=0,
#             edgecolors=c_score)
# ax_s.view_init(elev=22.5, azim=36)
# ax_s.set_xlabel(r'$lnb_1$')
# ax_s.set_ylabel(r'$lnr_b$')
# ax_s.set_zlabel(r'$f$')

# plt.show()

# # plot score in moment space
# # plt.close('all')

# fig = plt.figure(figsize=(6, 6))
# ax_s = fig.add_subplot(projection='3d')

# # connect points with the same ra
# ax_s.scatter(stats[index_p_all,0], stats[index_p_all,1], stats[index_p_all,2], 
#             'o',
#             s=100,
#             c = c_score,
#             alpha=1,
#             lw=0,
#             edgecolors=c_score)
# ax_s.view_init(elev=20, azim=36)
# ax_s.set_xlabel(r'$\mu$')
# ax_s.set_ylabel(r'$\sigma^{2}$')
# ax_s.set_zlabel(r'$\gamma$')

# plt.show()

#%% plot scree
fig = plt.figure(figsize=(6, 6))
ax_var = fig.add_subplot()
ax_var.plot(np.arange(len(S)),S)

ax_var.set_yscale('log')
ax_var.set_xlabel('rank')
ax_var.set_ylabel(r'$\Sigma$')

