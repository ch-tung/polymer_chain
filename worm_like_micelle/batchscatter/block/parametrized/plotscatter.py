# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 23:23:03 2021
plot scattering function
@author: CHTUNG
"""

import numpy as np
import numpy.matlib
from scipy.io import loadmat
# from scipy.io import savemat
# import matplotlib
import matplotlib.pyplot as plt
# from scipy import interpolate

class Plotscatter():
    def __init__(self):
        self.qq = 2*np.pi/(np.logspace(1,5,64))
        self.n_plot = 80
        self.n_color = 4
        
    def load(self,filename):
        scatter_dict = loadmat(filename)
        self.S_q  = scatter_dict['S_q']
        self.p = scatter_dict['p']
        
    def plot_setup(self):
        fig = plt.figure(figsize=(5, 5),dpi=192)
        self.ax = fig.add_subplot()
        x = np.linspace(0.0, 1.0, self.n_color)
        self.color = plt.get_cmap('viridis')(x)
        
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.set_xlim([5e-4, 2e-1])
        self.ax.set_ylim([0.5e-2, 2e0])
        self.ax.set_xlabel('$Q$')
        self.ax.set_ylabel('$S(Q)$')
        self.ax.grid(True,which='major')
        
        for i in range(16):
            self.ax.plot((10**(2*i+2)*np.array([1e-4, 1e1]))**-(1/4),np.array([1e-4, 1e1]),
                '--',color='#C0C0C0',linewidth=0.5)
            self.ax.plot((10**(i+1)*np.array([1e-4, 1e1]))**-(1/2),np.array([1e-4, 1e1]),
                    '--',color='#C0C0C0',linewidth=0.5)
            self.ax.plot((10**(i+1)*np.array([1e-4, 1e1]))**-(0.588),np.array([1e-4, 1e1]),
                    ':',color='#C0C0C0',linewidth=0.5)
            self.ax.plot((10**(i+1)*np.array([1e-4, 1e1]))**-(1),np.array([1e-4, 1e1]),
                    '-.',color='#C0C0C0',linewidth=0.5)
            
    def plot(self,linestyle_str):
        for i in range(self.n_plot):
            color_i = self.color[i,:]
            self.ax.plot(self.qq,self.S_q[:,i],linestyle_str,color = color_i)
            
        plt.show()
        
Plot = Plotscatter()


#%% load
# load mat files
# grep shape
filename = 'scatter_chain_block_0.mat'
scatter_dict = loadmat(filename)
qq_max = 2
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']


n_mat = 11
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

# fix the deviation high q limit
import Sk
q = qq
L = 1000
b = L*2
S_q_rod = Sk.S_rod(q,L)

for i in range(S_q.shape[1]):
    S_q_i = S_q[:,i]
    S_q_i[S_q_i<S_q_rod] = S_q_rod[S_q_i<S_q_rod]
    S_q[:,i] = S_q_i

qq = qq[qq<qq_max]

Plot.qq = scatter_dict['qq'][0,:]

Plot.n_plot = 11
Plot.n_color = Plot.n_plot
Plot.plot_setup()
Plot.ax.set_xlim([5e-4, 2e-1])

# #%%
# filename = 'scatter_chain_sfr_woSA.mat'
# scatter_dict = loadmat(filename)
# S_q_sfr  = scatter_dict['S_q']
# qq_sfr  = scatter_dict['qq']
# a_sfr  = scatter_dict['a'].T


# Plot.ax.plot(qq_sfr.T,S_q_sfr[:,[0]],'--',color = [0,0,0])
# Plot.ax.plot(qq_sfr.T,S_q_sfr[:,[3]],'-.',color = [0,0,0])

#%%prepare input
# property to be presented 0
p_0 = p[0]
pc_0 = (p_0-np.nanmin(p_0))/(np.nanmax(p_0)-np.nanmin(p_0))

# property to be presented 1
p_1 = p[1]/1000
pc_1 = (p_1-np.nanmin(p_1))/(np.nanmax(p_1)-np.nanmin(p_1))

# property to be presented 2
p_2 = p[2]
pc_2 = (p_2-np.nanmin(p_2))/(np.nanmax(p_2)-np.nanmin(p_2))

# load stats
filename_stats = 'stats_block.mat'
filename_parameters = 'parameters_block.mat'
stats_dict = loadmat(filename_stats)
parameters_dict = loadmat(filename_parameters)

stats = stats_dict['statistics']
p_m1 = stats[:,0]
p_m2 = stats[:,1]
p_m3 = stats[:,2]

set_m1 = sorted(set(p_m1))
set_m2 = sorted(set(p_m2))
set_m3 = sorted(set(p_m3))

index_p_m1 = (p_m1 == set_m1[10])
index_p_m2 = (p_m2 == set_m2[10])
index_p_m3 = (p_m3 == set_m3[5])
index_p = index_p_m3&index_p_m2
# index_p = np.arange(len(p[1])) # all datapoints, bool

Plot.S_q = S_q[:,index_p]
Plot.p = p[:,index_p]
Plot.plot('-')




# scatter_dict = loadmat(filename)

# S_q  = scatter_dict['S_q']
# # qq  = scatter_dict['qq'].T
# a  = scatter_dict['a'].T

# qq = 2*np.pi/(np.logspace(1,5,64))

# fig = plt.figure(figsize=(5, 5),dpi=192)
# ax = fig.add_subplot()

# n_plot = 6
# x = np.linspace(0.0, 1.0, n_plot)
# color = plt.get_cmap('viridis')(x)

# for i in range(16):
#     ax.plot((10**(2*i+2)*np.array([1e-4, 1e1]))**-(1/4),np.array([1e-4, 1e1]),
#         '--',color='#C0C0C0',linewidth=0.5)
#     ax.plot((10**(i+1)*np.array([1e-4, 1e1]))**-(1/2),np.array([1e-4, 1e1]),
#             '--',color='#C0C0C0',linewidth=0.5)
#     ax.plot((10**(i+1)*np.array([1e-4, 1e1]))**-(0.588),np.array([1e-4, 1e1]),
#             ':',color='#C0C0C0',linewidth=0.5)
#     ax.plot((10**(i+1)*np.array([1e-4, 1e1]))**-(1),np.array([1e-4, 1e1]),
#             '-.',color='#C0C0C0',linewidth=0.5)

# for i in range(n_plot):
#     color_i = color[i,:]
    
#     ax.plot(qq,S_q[:,i],'-',color = color_i)
    
# filename = 'scatter_chain_prstnc_woSA.mat'
# scatter_dict = loadmat(filename)

# for i in range(n_plot):
#     color_i = color[i,:]
    
#     ax.plot(qq,S_q[:,i],'--',color = color_i)
    
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim([np.sqrt(np.min(qq)*np.max(qq))*10**-2,np.sqrt(np.min(qq)*np.max(qq))*10**2])
# ax.set_ylim([2e-4, 2e0])
# ax.set_xlabel('$Q$')
# ax.set_ylabel('$S(Q)$')
# ax.grid(True,which='major')
# plt.show()