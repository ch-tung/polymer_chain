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

filename = 'scatter_chain.mat'
scatter_dict = loadmat(filename)

S_q  = scatter_dict['S_q']
qq  = scatter_dict['qq'].T
a  = scatter_dict['a'].T

fig = plt.figure(figsize=(5, 5),dpi=192)
ax = fig.add_subplot()

n_plot = 7
x = np.linspace(0.0, 1.0, n_plot)
color = plt.get_cmap('viridis')(x)

for i in range(5):
    ax.plot((10**(i+3)*np.array([1e-3, 1e1]))**-(1/2),np.array([1e-3, 1e1]),
            '--',color='#C0C0C0',linewidth=0.5)
    ax.plot((10**(i+3)*np.array([1e-3, 1e1]))**-(3/5),np.array([1e-3, 1e1]),
            ':',color='#C0C0C0',linewidth=0.5)

for i in range(n_plot):
    color_i = color[i,:]
    
    ax.plot(qq,S_q[:,i],'-',color = color_i)
    
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([np.sqrt(np.min(qq)*np.max(qq))*10**-1.5,np.sqrt(np.min(qq)*np.max(qq))*10**1.5])
ax.set_ylim([1e-3, 2e0])
ax.set_xlabel('$Q$')
ax.set_ylabel('$S(Q)$')
ax.grid(True,which='major')
plt.show()