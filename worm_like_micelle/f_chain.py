# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 12:27:38 2021
chain
@author: CHTUNG
"""
import numpy as np
import f_rotation
rotation = f_rotation.rotation
rotation_dihedral = f_rotation.rotation_dihedral
# rotation_stretched = f_rotation.rotation_stretched

def chain_Rayleigh(N, a, lambda_seg, unit_C, apply_SA=1, d_exc=1):
    d2_exc = d_exc**2
    i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg)) 
    # Check for sphere overlap was done for points 
    # separated by more than pi*d_exc/2 along the contour
       
    n = np.zeros((3,N))
    l = np.zeros((3,N))
    lc = np.zeros((3,N))
    #B = np.zeros((3,3))
    #C = np.zeros((3,3))
    #D = np.zeros((3,3))
    R = np.zeros((3,3))
    O = np.zeros((3,3,N))
    
    abort = 1
    while abort==1:
        abort = 0
        for i in range(N):
            if i==0:
                n[:,i] = [1,0,0]
                l[:,i] = n[:,i]
                #B = np.eye(3)
                #C = np.eye(3)
                #D = np.eye(3)
                R = np.eye(3)
                O[:,:,i] = R
            else:
                R = rotation(O[:,:,i-1],a)
                
                O[:,:,i] = R@O[:,:,i-1]
                # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                n[:,i] = O[:,0,i].reshape((3))
                # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<i_diameter:
                    continue
                
                #%% check self avoiding
                if apply_SA:
                    SA = 0
                    
                    n_retry = -1
                    while (SA == 0) & (n_retry < 100):
                        n_retry += 1
                        
                        # if n_retry > 100:
                        #     abort = 1
                        #     print('abort')
                        #     break
                            
                        d2_uv_min = np.min(np.sum((l[:,:i-i_diameter+1].T-l[:,i].T)**2,axis=1))
                        # d1_uv_min = np.min(np.max(np.abs(l[:,:i-1].T-l[:,i].T),axis=1))
                        # print(d1_uv_min)
                        
                        if d2_uv_min<d2_exc:
                        # if d1_uv_min<d_exc:
                            print('retry ({:d})'.format(n_retry+1))
                            # n_retry+=1
                            R = rotation(O[:,:,i-1],a)
                
                            O[:,:,i] = R@O[:,:,i-1]
                            # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                            n[:,i] = O[:,1,i].reshape((3))
                            # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                            l[:,i] = l[:,i-1] + n[:,i]
                        else:
                            if n_retry!=0:
                                print('retry (end)')
                            break
                        
                    if n_retry >= 100:
                        abort = 1
                        print('abort')
                        break
        
    lc = l*lambda_seg

    #%% map unimer
    #C
    nC = unit_C.shape[1]
    m_backbone_C = np.zeros((3,nC,N))
    for j in range(N):
        for k in range(nC):
            m_backbone_C[:,k,j] = O[:,:,j]@unit_C[:,k] + lc[:,j] + np.array([0,0,0])
    
    Cc = np.reshape(m_backbone_C,(3,N*nC))
    
    # print(n_retry)
    return lc, Cc, O, n

#%%
def chain_Rayleigh_woSA(N, a, lambda_seg, unit_C, apply_SA=0, d_exc=1):
    d2_exc = d_exc**2
    i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg)) 
    # Check for sphere overlap was done for points 
    # separated by more than pi*d_exc/2 along the contour
       
    n = np.zeros((3,N))
    l = np.zeros((3,N))
    lc = np.zeros((3,N))
    #B = np.zeros((3,3))
    #C = np.zeros((3,3))
    #D = np.zeros((3,3))
    R = np.zeros((3,3))
    O = np.zeros((3,3,N))
    
    abort = 1
    while abort==1:
        abort = 0
        for i in range(N):
            if i==0:
                n[:,i] = [1,0,0]
                l[:,i] = n[:,i]
                #B = np.eye(3)
                #C = np.eye(3)
                #D = np.eye(3)
                R = np.eye(3)
                O[:,:,i] = R
            else:
                R = rotation(O[:,:,i-1],a)
                
                O[:,:,i] = R@O[:,:,i-1]
                # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                n[:,i] = O[:,0,i].reshape((3))
                # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<i_diameter:
                    continue
        
    lc = l*lambda_seg

    #%% map unimer
    #C
    nC = unit_C.shape[1]
    m_backbone_C = np.zeros((3,nC,N))
    for j in range(N):
        for k in range(nC):
            m_backbone_C[:,k,j] = O[:,:,j]@unit_C[:,k] + lc[:,j] + np.array([0,0,0])
    
    Cc = np.reshape(m_backbone_C,(3,N*nC))
    
    # print(n_retry)
    return lc, Cc, O, n

#%%
def chain_fix_val_free_rot(N, a, lambda_seg, unit_C, apply_SA=1, d_exc=1):
    d2_exc = d_exc**2
    i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg)) 
    # Check for sphere overlap was done for points 
    # separated by more than pi*d_exc/2 along the contour
       
    n = np.zeros((3,N))
    l = np.zeros((3,N))
    lc = np.zeros((3,N))
    #B = np.zeros((3,3))
    #C = np.zeros((3,3))
    #D = np.zeros((3,3))
    R = np.zeros((3,3))
    O = np.zeros((3,3,N))
    
    abort = 1
    while abort==1:
        abort = 0
        for i in range(N):
            if i==0:
                n[:,i] = [1,0,0]
                l[:,i] = n[:,i]
                #B = np.eye(3)
                #C = np.eye(3)
                #D = np.eye(3)
                R = np.eye(3)
                O[:,:,i] = R
            else:
                R = rotation_dihedral(O[:,:,i-1],a)
                
                O[:,:,i] = R@O[:,:,i-1]
                # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                n[:,i] = O[:,0,i].reshape((3))
                # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<i_diameter:
                    continue
                
                #%% check self avoiding
                if apply_SA:
                    SA = 0
                    
                    n_retry = -1
                    while (SA == 0) & (n_retry < 100):
                        n_retry += 1
                        
                        # if n_retry > 100:
                        #     abort = 1
                        #     print('abort')
                        #     break
                            
                        d2_uv_min = np.min(np.sum((l[:,:i-i_diameter+1].T-l[:,i].T)**2,axis=1))
                        # d1_uv_min = np.min(np.max(np.abs(l[:,:i-1].T-l[:,i].T),axis=1))
                        # print(d1_uv_min)
                        
                        if d2_uv_min<d2_exc:
                        # if d1_uv_min<d_exc:
                            print('retry ({:d})'.format(n_retry+1))
                            # n_retry+=1
                            R = rotation(O[:,:,i-1],a)
                
                            O[:,:,i] = R@O[:,:,i-1]
                            # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                            n[:,i] = O[:,1,i].reshape((3))
                            # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                            l[:,i] = l[:,i-1] + n[:,i]
                        else:
                            if n_retry!=0:
                                print('retry (end)')
                            break
                        
                    if n_retry >= 100:
                        abort = 1
                        print('abort')
                        break
        
    lc = l*lambda_seg
    
    #%% map unimer
    #C
    nC = unit_C.shape[1]
    m_backbone_C = np.zeros((3,nC,N))
    for j in range(N):
        for k in range(nC):
            m_backbone_C[:,k,j] = O[:,:,j]@unit_C[:,k] + lc[:,j] + np.array([0,0,0])
    
    Cc = np.reshape(m_backbone_C,(3,N*nC))
    
    # print(n_retry)
    return lc, Cc, O, n

#%%
def chain_fix_val_free_rot_woSA(N, a, lambda_seg, unit_C, apply_SA=0, d_exc=1):
    d2_exc = d_exc**2
    i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg)) 
    # Check for sphere overlap was done for points 
    # separated by more than pi*d_exc/2 along the contour
       
    n = np.zeros((3,N))
    l = np.zeros((3,N))
    lc = np.zeros((3,N))
    #B = np.zeros((3,3))
    #C = np.zeros((3,3))
    #D = np.zeros((3,3))
    R = np.zeros((3,3))
    O = np.zeros((3,3,N))
    
    abort = 1
    while abort==1:
        abort = 0
        for i in range(N):
            if i==0:
                n[:,i] = [1,0,0]
                l[:,i] = n[:,i]
                #B = np.eye(3)
                #C = np.eye(3)
                #D = np.eye(3)
                R = np.eye(3)
                O[:,:,i] = R
            else:
                R = rotation_dihedral(O[:,:,i-1],a)
                
                O[:,:,i] = R@O[:,:,i-1]
                # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                n[:,i] = O[:,0,i].reshape((3))
                # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<i_diameter:
                    continue
        
    lc = l*lambda_seg
    
    #%% map unimer
    #C
    nC = unit_C.shape[1]
    m_backbone_C = np.zeros((3,nC,N))
    for j in range(N):
        for k in range(nC):
            m_backbone_C[:,k,j] = O[:,:,j]@unit_C[:,k] + lc[:,j] + np.array([0,0,0])
    
    Cc = np.reshape(m_backbone_C,(3,N*nC))
    
    # print(n_retry)
    return lc, Cc, O, n

#%% 
import random
def chain_grid(N, kappa, epsilon, lambda_seg, apply_SA=1, d_exc=1):
    d2_exc = d_exc**2
    i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg))
    
    # coordinate 
    # 6 orientations following the arrangement of standard dice
    # 1/6 for +x/-x; 2/5 for +y/-y; 3/4 for +z/-z
    r_n = np.array([[1,0,0],[0,1,0],[0,0,1],[0,0,-1],[0,-1,0],[-1,0,0]])
    r_opp = np.array([5,4,3,2,1,0])
    
    cos_ij = r_n@r_n.T
    sin_ij2 = 1-cos_ij**2
    
    n = np.zeros((3,N))
    n_i = np.zeros((N)).astype(int)
    l = np.zeros((3,N))
    
    # energy
    Z = np.zeros((6,6))
    for i in range(6):
        E_phi = kappa*(sin_ij2[i,:])
        E_x = -epsilon*(cos_ij[0,:])
        
        E = E_phi + E_x
        z_i = np.exp(-E)
        z_i[r_opp[i]] = 0     
        
        z_i = z_i/np.sum(z_i)
        Z[:,i] = z_i
    
    abort = 1
    while abort==1:
        abort = 0
        for i in range(N):
            if i==0:
                n_i[i] = 0
                n[:,i] = r_n[n_i[i],:]
                l[:,i] = n[:,i]
                
            else:
                z = Z[:,n_i[i-1]]
                n_i[i] = random.choices(np.arange(6),weights=z)[0]
                            
                n[:,i] = r_n[n_i[i],:]
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<i_diameter:
                        continue
                
                #%% check self avoiding
                if apply_SA:
                    SA = 0
                    
                    n_retry = -1
                    while (SA == 0) & (n_retry < 100):
                        n_retry += 1
                        
                        # if n_retry > 100:
                        #     abort = 1
                        #     print('abort')
                        #     break
                            
                        d2_uv_min = np.min(np.sum((l[:,:i-i_diameter+1].T-l[:,i].T)**2,axis=1))
                        # d1_uv_min = np.min(np.max(np.abs(l[:,:i-1].T-l[:,i].T),axis=1))
                        # print(d1_uv_min)
                        
                        if d2_uv_min<d2_exc:
                        # if d1_uv_min<d_exc:
                            print('retry ({:d})'.format(n_retry+1))
                            # n_retry+=1
                            z = Z[:,n_i[i-1]]
                            n_i[i] = random.choices(np.arange(6),weights=z)[0]
                                        
                            n[:,i] = r_n[n_i[i],:]
                            l[:,i] = l[:,i-1] + n[:,i]
                            
                        else:
                            if n_retry!=0:
                                print('retry (end)')
                            break
                        
                    if n_retry >= 100:
                        abort = 1
                        print('abort')
                        break
                    
    lc = l*lambda_seg
    Cc = lc
    
    return lc, Cc, n

def chain_grid_woSA(N, kappa, epsilon, lambda_seg, apply_SA=1, d_exc=0):
    # d2_exc = d_exc**2
    i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg))
    
    # coordinate 
    # 6 orientations following the arrangement of standard dice
    # 1/6 for +x/-x; 2/5 for +y/-y; 3/4 for +z/-z
    r_n = np.array([[1,0,0],[0,1,0],[0,0,1],[0,0,-1],[0,-1,0],[-1,0,0]])
    r_opp = np.array([5,4,3,2,1,0])
    
    cos_ij = r_n@r_n.T
    sin_ij2 = 1-cos_ij**2
    
    n = np.zeros((3,N))
    n_i = np.zeros((N)).astype(int)
    l = np.zeros((3,N))
    
    # energy
    Z = np.zeros((6,6))
    for i in range(6):
        E_phi = kappa*(sin_ij2[i,:])
        E_x = -epsilon*(cos_ij[0,:])
        
        E = E_phi + E_x
        z_i = np.exp(-E)
        z_i[r_opp[i]] = 0     
        
        z_i = z_i/np.sum(z_i)
        Z[:,i] = z_i
    
    abort = 1
    while abort==1:
        abort = 0
        for i in range(N):
            if i==0:
                n_i[i] = 0
                n[:,i] = r_n[n_i[i],:]
                l[:,i] = n[:,i]
                
            else:
                z = Z[:,n_i[i-1]]
                n_i[i] = random.choices(np.arange(6),weights=z)[0]
                            
                n[:,i] = r_n[n_i[i],:]
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<i_diameter:
                        continue
                    
    lc = l*lambda_seg
    Cc = lc
    
    return lc, Cc, n
#%%
# def chain_stretched(N, a, lambda_seg, unit_C, apply_SA=1, d_exc=1):
#     d2_exc = d_exc**2
#     i_diameter = int(np.ceil(np.pi/2*d_exc/lambda_seg)) 
#     # Check for sphere overlap was done for points 
#     # separated by more than pi*d_exc/2 along the contour
       
#     n = np.zeros((3,N))
#     l = np.zeros((3,N))
#     lc = np.zeros((3,N))
#     #B = np.zeros((3,3))
#     #C = np.zeros((3,3))
#     #D = np.zeros((3,3))
#     R = np.zeros((3,3))
#     O = np.zeros((3,3,N))
    
#     abort = 1
#     while abort==1:
#         abort = 0
#         for i in range(N):
#             if i==0:
#                 n[:,i] = [1,0,0]
#                 l[:,i] = n[:,i]
#                 #B = np.eye(3)
#                 #C = np.eye(3)
#                 #D = np.eye(3)
#                 R = np.eye(3)
#                 O[:,:,i] = R
#             else:
#                 R = rotation(O[:,:,i-1],a)
                
#                 O[:,:,i] = R@O[:,:,i-1]
#                 # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
#                 n[:,i] = O[:,0,i].reshape((3))
#                 # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
#                 l[:,i] = l[:,i-1] + n[:,i]
                
#                 if i<i_diameter:
#                     continue
                
#                 #%% check self avoiding
#                 if apply_SA:
#                     SA = 0
                    
#                     n_retry = -1
#                     while (SA == 0) & (n_retry < 100):
#                         n_retry += 1
                        
#                         # if n_retry > 100:
#                         #     abort = 1
#                         #     print('abort')
#                         #     break
                            
#                         d2_uv_min = np.min(np.sum((l[:,:i-i_diameter+1].T-l[:,i].T)**2,axis=1))
#                         # d1_uv_min = np.min(np.max(np.abs(l[:,:i-1].T-l[:,i].T),axis=1))
#                         # print(d1_uv_min)
                        
#                         if d2_uv_min<d2_exc:
#                         # if d1_uv_min<d_exc:
#                             print('retry ({:d})'.format(n_retry+1))
#                             # n_retry+=1
#                             R = rotation(O[:,:,i-1],a)
                
#                             O[:,:,i] = R@O[:,:,i-1]
#                             # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
#                             n[:,i] = O[:,1,i].reshape((3))
#                             # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
#                             l[:,i] = l[:,i-1] + n[:,i]
#                         else:
#                             if n_retry!=0:
#                                 print('retry (end)')
#                             break
                        
#                     if n_retry >= 100:
#                         abort = 1
#                         print('abort')
#                         break
        
#     lc = l*lambda_seg

#     #%% map unimer
#     #C
#     nC = unit_C.shape[1]
#     m_backbone_C = np.zeros((3,nC,N))
#     for j in range(N):
#         for k in range(nC):
#             m_backbone_C[:,k,j] = O[:,:,j]@unit_C[:,k] + lc[:,j] + np.array([0,0,0])
    
#     Cc = np.reshape(m_backbone_C,(3,N*nC))
    
#     # print(n_retry)
#     return lc, Cc, O, n