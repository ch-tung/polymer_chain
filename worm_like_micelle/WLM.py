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
def rotation(O,a):
    # quaternion
    phi_q = 2*(np.random.rand(1)-0.5)*np.pi
    delta = 1e-3
    theta_q = np.sqrt(-np.log(1-np.random.rand(1)*(1-delta))/a);
    if theta_q>np.pi/3*2:
        theta_q = np.array([np.pi/3*2])
    #print(theta_q)
    
    vq = O[:,1]*np.cos(phi_q) + O[:,2]*np.sin(phi_q)
    vq = vq/np.sqrt(np.sum(vq**2))
    qr = np.cos(theta_q/2);
    qi = vq[0]*np.sin(theta_q/2);
    qj = vq[1]*np.sin(theta_q/2);
    qk = vq[2]*np.sin(theta_q/2);
    nq = np.sqrt(qr**2 + qi**2 + qj**2 + qk**2)
    qr = qr/nq
    qi = qi/nq
    qj = qj/nq
    qk = qk/nq
    
    Rq = np.array([[1-2*(qj**2+qk**2), 2*(qi*qj+qk*qr), 2*(qi*qk-qj*qr)],
                    [2*(qi*qj-qk*qr), 1-2*(qi**2+qk**2), 2*(qj*qk+qi*qr)],
                    [2*(qi*qk+qj*qr), 2*(qj*qk-qi*qr), 1-2*(qj**2+qk**2)]])
    #print(Rq)
    
    R = Rq[:,:,0]
    return R
    
def chain_Rayleigh(DP, a, lambda_seg, unit_C, apply_SA=True, d_exc=0.5):
       
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
            # # quaternion
            # phi_q = 2*(np.random.rand(1)-0.5)*np.pi
            # delta = 1e-4
            # theta_q = np.sqrt(-np.log(1-np.random.rand(1)*(1-delta))/a);
            # if theta_q>np.pi:
            #     theta_q = np.array([np.pi-delta])
            # #print(theta_q)
            
            # vq = O[:,1,i-1]*np.cos(phi_q) + O[:,2,i-1]*np.sin(phi_q)
            # qr = np.cos(theta_q/2);
            # qi = vq[0]*np.sin(theta_q/2);
            # qj = vq[1]*np.sin(theta_q/2);
            # qk = vq[2]*np.sin(theta_q/2);
            
            # Rq = np.array([[1-2*(qj**2+qk**2), 2*(qi*qj+qk*qr), 2*(qi*qk-qj*qr)],
            #                 [2*(qi*qj-qk*qr), 1-2*(qi**2+qk**2), 2*(qj*qk+qi*qr)],
            #                 [2*(qi*qk+qj*qr), 2*(qj*qk-qi*qr), 1-2*(qj**2+qk**2)]])
            # #print(Rq)
            
            R = rotation(O[:,:,i-1],a)
            
            O[:,:,i] = R@O[:,:,i-1]
            O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
            n[:,i] = O[:,1,i].reshape((3))
            # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
            l[:,i] = l[:,i-1] + n[:,i]
            
            #%% check self avoiding
            if apply_SA:
                for u in range(i-1):
                    d_uv = np.sqrt(np.sum((l[:,i-u-1]-l[:,i])**2))
                    if d_uv<d_exc:
                        #print('retry')
                        R = rotation(O[:,:,i-1],a)
            
                        O[:,:,i] = R@O[:,:,i-1]
                        O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                        n[:,i] = O[:,1,i].reshape((3))
                        # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                        l[:,i] = l[:,i-1] + n[:,i]
                    else:
                        break
            
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
DP_backbone = 10000

# Chain stiffness
a_backbone = 1e2

# Unit persistence
lambda_backbone = 1

# call 'chain' function
lc_backbone, Cc_backbone, O_backbone, n_backbone = chain_Rayleigh(DP_backbone,a_backbone,lambda_backbone,unit_C)

tEnd = time.time()
print("It cost %f sec" % (tEnd - tStart))
#%% plot
fig = plt.figure(figsize=(6, 6),dpi=192)
ax = fig.add_subplot(projection='3d')

ax.plot(Cc_backbone[0,:],Cc_backbone[1,:],Cc_backbone[2,:], 
        '-', color='#D00000', linewidth=2, alpha = 0.75)
# ax.plot(Cc_backbone[0,:],Cc_backbone[1,:],Cc_backbone[2,:], 
#         'o', markeredgecolor='#800000', markerfacecolor='#D00000')

ax.plot(Cc_backbone[0,0],Cc_backbone[1,0],Cc_backbone[2,0], 
            'o', markeredgecolor='#800000', markerfacecolor='#D00000')
ax.plot(Cc_backbone[0,-1],Cc_backbone[1,-1],Cc_backbone[2,-1], 
            'o', markeredgecolor='#800000', markerfacecolor='#D00000')

CM = np.mean(Cc_backbone,axis=1)
CT = np.array([np.max(Cc_backbone[0,:])+np.min(Cc_backbone[0,:]),
               np.max(Cc_backbone[1,:])+np.min(Cc_backbone[1,:]),
               np.max(Cc_backbone[2,:])+np.min(Cc_backbone[2,:])])/2
d_box = np.max([np.max(Cc_backbone[0,:])-np.min(Cc_backbone[0,:]),
                np.max(Cc_backbone[1,:])-np.min(Cc_backbone[1,:]),
                np.max(Cc_backbone[2,:])-np.min(Cc_backbone[2,:])])

#ax.axis('off')
ax.set_xlim([CT[0]-d_box/2, CT[0]+d_box/2])
ax.set_ylim([CT[1]-d_box/2, CT[1]+d_box/2])
ax.set_zlim([CT[2]-d_box/2, CT[2]+d_box/2])
ax.set_box_aspect([1,1,1])
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.set_zticklabels([])

plt.show()

#%%
# fig2 = plt.figure(figsize=(6, 6),dpi=192)
# ax2 = fig2.add_subplot()

# ax2.plot(np.arange(DP_backbone), np.sum(n_backbone**2,axis=0),'-')
# ax2.set_yscale('log')
# ax2.set_ylim([0.5, 2])
