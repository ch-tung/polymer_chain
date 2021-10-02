# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 22:23:28 2021
Generate WLM chain trajectories
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
from scipy.io import loadmat
from scipy.io import savemat
import matplotlib.pyplot as plt

import time
tStart = time.time()

#%% define functions

def chain_Rayleigh(DP,a,lambda_seg,unit_C):
    n = np.zeros((3,DP))
    l = np.zeros((3,DP))
    lc = np.zeros((3,DP))
    #B = np.zeros((3,3))
    #C = np.zeros((3,3))
    #D = np.zeros((3,3))
    R = np.zeros((3,3))
    O = np.zeros((3,3,DP))
    
    for i in range(DP):
        if i==0:
            n[:,i] = [1,0,0]
            l[:,i] = n[:,i]
            #B = np.eye(3)
            #C = np.eye(3)
            #D = np.eye(3)
            R = np.eye(3)
            O[:,:,i] = R
        else:
            # quaternion
            phi_q = 2*(np.random.rand(1)-0.5)*np.pi
            theta_q = np.sqrt(-np.log(1-np.random.rand(1))/a);
            
            vq = O[:,1,i-1]*np.cos(phi_q) + O[:,2,i-1]*np.sin(phi_q)
            qr = np.cos(theta_q/2);
            qi = vq[0]*np.sin(theta_q/2);
            qj = vq[1]*np.sin(theta_q/2);
            qk = vq[2]*np.sin(theta_q/2);
            
            Rq = np.array([[1-2*(qj**2+qk**2), 2*(qi*qj+qk*qr), 2*(qi*qk-qj*qr)],
                            [2*(qi*qj-qk*qr), 1-2*(qi**2+qk**2), 2*(qj*qk+qi*qr)],
                            [2*(qi*qk+qj*qr), 2*(qj*qk-qi*qr), 1-2*(qj**2+qk**2)]])
            
            R = Rq[:,:,0]
            
            O[:,:,i] = R@O[:,:,i-1]
            n[:,i] = R@n[:,i-1]
            l[:,i] = l[:,i-1] + n[:,i]
            
    lc = l*lambda_seg

    #%% map unimer
    #C
    nC = unit_C.shape[1]
    m_backbone_C = np.zeros((3,nC,DP))
    for j in range(DP):
        for k in range(nC):
            m_backbone_C[:,k,j] = O[:,:,j]@unit_C[:,k] + lc[:,j] + np.array([0,0,0])
    
    Cc = np.reshape(m_backbone_C,(3,DP*nC))
    
    return lc, Cc, O, n

#%% backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
DP_backbone = 100

# Chain stiffness
a_backbone = 10

# Unit persistence
lambda_backbone = 1

# call 'chain' function
lc_backbone, Cc_backbone, O_backbone, n_backbone = chain_Rayleigh(DP_backbone,a_backbone,lambda_backbone,unit_C)

#%% plot
fig = plt.figure(figsize=(6, 6),dpi=192)
ax = fig.add_subplot(projection='3d')

ax.plot(Cc_backbone[0,:],Cc_backbone[1,:],Cc_backbone[2,:], 
           '-', color='#303030', linewidth=2, markersize=12)
ax.plot(Cc_backbone[0,:],Cc_backbone[1,:],Cc_backbone[2,:], 
           'o', markeredgecolor='#800000', markerfacecolor='#D00000')

ax.axis('off')

plt.show()