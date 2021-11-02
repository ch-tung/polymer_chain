# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 23:23:03 2021
plot scattering function
@author: CHTUNG
"""

import numpy as np
import numpy.matlib
from scipy.io import loadmat
from scipy.io import savemat
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate

class Plotscatter():
    def __init__(self):
        self.qq = 2*np.pi/(np.logspace(1,5,64))
        self.n_plot = 6
        
    def load(self,filename):
        scatter_dict = loadmat(filename)
        self.S_q  = scatter_dict['S_q']
        self.a = a  = scatter_dict['a'].T
        
    def plot_setup(self):
        fig = plt.figure(figsize=(5, 5),dpi=192)
        self.ax = fig.add_subplot()
        n_plot = 6
        x = np.linspace(0.0, 1.0, n_plot)
        self.color = plt.get_cmap('viridis')(x)
        
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.set_xlim([np.sqrt(np.min(self.qq)*np.max(self.qq))*10**-2,np.sqrt(np.min(self.qq)*np.max(self.qq))*10**2])
        self.ax.set_ylim([2e-4, 2e0])
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
Plot.plot_setup()

filename = 'scatter_chain_prstnc_woSA.mat'
Plot.load(filename)
Plot.plot('--')

filename = 'scatter_chain_prstnc.mat'
Plot.load(filename)
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