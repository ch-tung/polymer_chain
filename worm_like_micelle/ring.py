# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 20:54:38 2021
Ring
@author: CHTUNG
"""
import sys
sys.modules[__name__].__dict__.clear()

import numpy as np
import numpy.matlib
#from scipy.io import loadmat
#from scipy.io import savemat
import matplotlib.pyplot as plt
from WLM import WLChain
import time

# parameters
N = 1024
n_harmonics = 10

# c_ring = np.zeros((3,N+1))

# for i in range(3):
#     phi_i = 2*np.pi*np.random.rand(1)
    
#     weight = 1/(np.arange(n_harmonics)+1)**1.5
#     weight = weight/np.sqrt(np.sum(weight**2))
#     coeff_c_i = np.random.rand(n_harmonics)*weight
#     coeff_s_i = np.random.rand(n_harmonics)*weight
    
#     theta = np.arange(N+1)/N*2*np.pi
    
#     harmonics_c_i = np.cos(np.outer(theta,(np.arange(n_harmonics)+1)) + phi_i)*coeff_c_i
#     harmonics_s_i = np.sin(np.outer(theta,(np.arange(n_harmonics)+1)) + phi_i)*coeff_s_i
    
#     harmonics_i = harmonics_c_i + harmonics_s_i
    
#     c_ring[i,:] = np.sum(harmonics_i,axis=1)

chain01 = WLChain(N=1024)
plt.close('all')
for i in range(20):
    tStart = time.time()
    #chain01.apply_SA = 0
    chain01.ring(n_harmonics=10)
    tEnd = time.time()
    
    filename_ring = './figures/ring/ring_test_{:d}.png'.format(i+1)
    chain01.plot(filename=filename_ring, show_axes=0, save=1, end=0)
#%% plot


# show_axes=0
# fig = plt.figure(figsize=(6, 6),dpi=192)
# ax = fig.add_subplot(projection='3d')

# ax.plot(c_ring[0,:],c_ring[1,:],c_ring[2,:], 
#         '-', color='#D00000', linewidth=2, alpha = 0.75)
# # ax.plot(c_ring[0,:],c_ring[1,:],c_ring[2,:], 
# #         'o', markeredgecolor='#800000', markerfacecolor='#D00000')

# # ax.plot(c_ring[0,0],c_ring[1,0],c_ring[2,0], 
# #             'o', markeredgecolor='#800000', markerfacecolor='#D00000')
# # ax.plot(c_ring[0,-1],c_ring[1,-1],c_ring[2,-1], 
# #             'o', markeredgecolor='#800000', markerfacecolor='#D00000')

# #CM = np.mean(Cc_backbone,axis=1)
# CT = np.array([np.max(c_ring[0,:])+np.min(c_ring[0,:]),
#                 np.max(c_ring[1,:])+np.min(c_ring[1,:]),
#                 np.max(c_ring[2,:])+np.min(c_ring[2,:])])/2
# d_box = np.max([np.max(c_ring[0,:])-np.min(c_ring[0,:]),
#                 np.max(c_ring[1,:])-np.min(c_ring[1,:]),
#                 np.max(c_ring[2,:])-np.min(c_ring[2,:])])

# if show_axes==0:
#     ax.set_xticklabels([])
#     ax.set_yticklabels([])
#     ax.set_zticklabels([])
#     #ax.axis('off')

# ax.set_xlim([CT[0]-d_box/2, CT[0]+d_box/2])
# ax.set_ylim([CT[1]-d_box/2, CT[1]+d_box/2])
# ax.set_zlim([CT[2]-d_box/2, CT[2]+d_box/2])
# ax.set_box_aspect([1,1,1])

# plt.show()