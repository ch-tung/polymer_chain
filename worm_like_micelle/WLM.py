# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 22:23:28 2021
Generate WLM chain trajectories
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
#from scipy.io import loadmat
#from scipy.io import savemat
import matplotlib.pyplot as plt

#import time

#%% define functions
def rotation(O,a):
    # quaternion
    phi_q = 2*(np.random.rand(1)-0.5)*np.pi
    delta = 1e-3
    theta_q = np.sqrt(-np.log(1-np.random.rand(1)*(1-delta))/a)/2;
    # if theta_q>np.pi/3*2:
    #     theta_q = np.array([np.pi/3*2])/2
    #print(theta_q)
    sin_theta_q = np.sin(theta_q)
    
    vq = O[:,1]*np.cos(phi_q) + O[:,2]*np.sin(phi_q)
    # vq = vq/np.sqrt(np.sum(vq**2))
    qr = np.cos(theta_q);
    qi = vq[0]*sin_theta_q;
    qj = vq[1]*sin_theta_q;
    qk = vq[2]*sin_theta_q;
    # nq = np.sqrt(qr**2 + qi**2 + qj**2 + qk**2)
    # qr = qr/nq
    # qi = qi/nq
    # qj = qj/nq
    # qk = qk/nq
    
    qij = qi*qj
    qjk = qj*qk
    qik = qi*qk
    qir = qi*qr
    qjr = qj*qr
    qkr = qk*qr
    qii = qi*qi
    qjj = qj*qj
    qkk = qk*qk
    
    Rq = np.array([[1-2*(qjj+qkk), 2*(qij+qkr), 2*(qik-qjr)],
                   [2*(qij-qkr), 1-2*(qii+qkk), 2*(qjk+qir)],
                   [2*(qik+qjr), 2*(qjk-qir), 1-2*(qii+qjj)]])
    
    R = Rq[:,:,0]
    
    # Re-orthogonalize
    Rx = R[:,0]
    Ry = R[:,1]
    err = np.dot(Rx,Ry)
    Rx_ort = Rx-(err/2)*Ry
    Ry_ort = Ry-(err/2)*Rx
    #Rz_ort = np.cross(Rx_ort,Ry_ort)
    Rx_new = 0.5*(3-np.dot(Rx_ort,Rx_ort))*Rx_ort
    Ry_new = 0.5*(3-np.dot(Ry_ort,Ry_ort))*Ry_ort
    Rz_new = np.cross(Rx_new,Ry_new)
    #Rz_new = 0.5*(3-np.dot(Rz_ort,Rz_ort))*Rz_ort
    R = np.array([Rx_new, Ry_new, Rz_new]).T
    
    return R
   
def chain_Rayleigh(N, a, lambda_seg, unit_C, apply_SA=1, d_exc=1):
    d2_exc = d_exc**2
       
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
                n[:,i] = O[:,1,i].reshape((3))
                # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                l[:,i] = l[:,i-1] + n[:,i]
                
                if i<2:
                    continue
                
                #%% check self avoiding
                if apply_SA:
                    SA = 0
                    
                    n_retry = -1
                    while SA == 0:
                        n_retry += 1
                        
                        if n_retry > 100:
                            abort = 1
                            print('abort')
                            break
                            
                        d2_uv_min = np.min(np.sum((l[:,:i-1].T-l[:,i].T)**2,axis=1))
                        
                        if d2_uv_min<d2_exc:                      
                            print('retry')
                            # n_retry+=1
                            R = rotation(O[:,:,i-1],a)
                
                            O[:,:,i] = R@O[:,:,i-1]
                            # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                            n[:,i] = O[:,1,i].reshape((3))
                            # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                            l[:,i] = l[:,i-1] + n[:,i]
                        else:
                            break
                        
                    if abort==1:
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

def ring_harmonic(N,n_harmonics):
    c_ring = np.zeros((3,N+1))
    for i in range(3):
        phi_i = 2*np.pi*np.random.rand(1)
        
        weight = np.exp(-(np.arange(n_harmonics)+1)**2/5)
        weight = weight/np.sqrt(np.sum(weight**2))
        coeff_c_i = np.random.rand(n_harmonics)*weight
        coeff_s_i = np.random.rand(n_harmonics)*weight
        
        theta = np.arange(N+1)/N*2*np.pi
        
        harmonics_c_i = np.cos(np.outer(theta,(np.arange(n_harmonics)+1)) + phi_i)*coeff_c_i
        harmonics_s_i = np.sin(np.outer(theta,(np.arange(n_harmonics)+1)) + phi_i)*coeff_s_i
        
        harmonics_i = harmonics_c_i + harmonics_s_i
        
        c_ring[i,:] = np.sum(harmonics_i,axis=1)
    return c_ring
#%% class: WLChain
class WLChain:
    """
    Modelling polymer chain based on the worm-like chain model and calculate their scattering function.
    """
    lc = []
    Cc = []
    O = []
    n = []
    l_contour = []
    l_end2end = []
    box = []
    apply_SA = []
    d_exc = []
    
    def __init__(self, N=1000, a=1e2, lmbda=1, unit_C=np.zeros((3,1))):
        self.N = N
        self.a = a
        self.lmbda = lmbda
        self.unit_C = unit_C
        self.apply_SA = 1
        self.d_exc = 1
        
    def chain(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_Rayleigh' function
        self.lc, self.Cc, self.O, self.n = chain_Rayleigh(self.N,self.a,self.lmbda,self.unit_C,
                                                          apply_SA=self.apply_SA,d_exc=self.d_exc)
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def ring(self,n_harmonics):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'ring_harmonics' function
        self.Cc = ring_harmonic(self.N,n_harmonics)
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
    
    def plot(self, filename=[], show_axes=1, save=0, end=1):
        """
        Plot polymer chain.
        
        Args:
            show_axes: boolean
            save: boolean
        """
        
        #plt.close('all')
        fig = plt.figure(figsize=(6, 6),dpi=192)
        ax = fig.add_subplot(projection='3d')
        
        ax.plot(self.Cc[0,:],self.Cc[1,:],self.Cc[2,:], 
                '-', color='#D00000', linewidth=2, alpha = 0.75)
        # ax.plot(self.Cc[0,:],self.Cc[1,:],self.Cc[2,:], 
        #         'o', markeredgecolor='#800000', markerfacecolor='#D00000')
        
        if end==1:
            ax.plot(self.Cc[0,0],self.Cc[1,0],self.Cc[2,0], 
                        'o', markeredgecolor='#800000', markerfacecolor='#D00000')
            ax.plot(self.Cc[0,-1],self.Cc[1,-1],self.Cc[2,-1], 
                        'o', markeredgecolor='#800000', markerfacecolor='#D00000')
        
        #CM = np.mean(Cc_backbone,axis=1)
        CT = np.array([np.max(self.Cc[0,:])+np.min(self.Cc[0,:]),
                       np.max(self.Cc[1,:])+np.min(self.Cc[1,:]),
                       np.max(self.Cc[2,:])+np.min(self.Cc[2,:])])/2
        d_box = np.max([np.max(self.Cc[0,:])-np.min(self.Cc[0,:]),
                        np.max(self.Cc[1,:])-np.min(self.Cc[1,:]),
                        np.max(self.Cc[2,:])-np.min(self.Cc[2,:])])
        
        if show_axes==0:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            #ax.axis('off')
        
        ax.set_xlim([CT[0]-d_box/2, CT[0]+d_box/2])
        ax.set_ylim([CT[1]-d_box/2, CT[1]+d_box/2])
        ax.set_zlim([CT[2]-d_box/2, CT[2]+d_box/2])
        ax.set_box_aspect([1,1,1])
        
        if save==1:
            plt.savefig(filename)

        plt.show()
        
    def scatter(self, n_grid=256, approx_1D=0):
        """
        Calculate scattering function.
        
        Args:
            n_grid: int
                number of grid points
            approx_1D: boolean
                1-D FFT for isotropic systems
        """
        
        N = self.N
        chain_box = self.box
    
        #box_size = np.max(chain_box[1,:]-chain_box[0,:],axis=0)
        box_size = N
        grid_size = (box_size)/n_grid
        Cc_relative = self.Cc.T-chain_box[0,:] # relative position of WL-chain in the box
        bead_coord = np.floor(Cc_relative/grid_size).astype('int')
        
        if approx_1D==1:
            # density in real space
            rho_rx = np.zeros(n_grid)
            rho_ry = np.zeros(n_grid)
            rho_rz = np.zeros(n_grid)
            
            for i in range(N):
                rho_rx[bead_coord[i,0]] += 1
                rho_ry[bead_coord[i,0]] += 1
                rho_rz[bead_coord[i,0]] += 1
            
            # FFT and calculate scattering function
            rho_qx = np.fft.fft(rho_rx)
            rho_qy = np.fft.fft(rho_ry)
            rho_qz = np.fft.fft(rho_rz)
            S_q_x = np.absolute(rho_qx)**2/N
            S_q_y = np.absolute(rho_qy)**2/N
            S_q_z = np.absolute(rho_qz)**2/N
            S_q_ave = (S_q_x + S_q_y + S_q_z)/3
            
            # radial average
            grid_coord = np.meshgrid(np.arange(n_grid))
            dq_grid = 2*np.pi/(box_size)
            q_grid = np.sqrt(grid_coord[0]**2)*dq_grid
            
            dq = dq_grid
            nq = int(np.floor(dq_grid/dq*n_grid/2))
            qq = (np.arange(nq)+0.5)*dq    
            index_q = np.floor(q_grid/dq) # q to the origin
            
            S_q = np.zeros(int(nq))
            
            for iq in range(int(nq)):
                S_q[iq] = np.sum(S_q_ave[index_q==iq])/N
                S_q[iq] = np.average(S_q_ave[index_q==iq])/N
            
        else:
            # density in real space
            rho_r = np.zeros((n_grid,n_grid,n_grid))
            
            for i in range(N):
                rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1
            
            # FFT and calculate scattering function
            rho_q = np.fft.fftn(rho_r)
            S_q_lmn = np.absolute(rho_q)**2/N
            
            # radial average
            grid_coord = np.meshgrid(np.arange(n_grid),np.arange(n_grid),np.arange(n_grid))
            dq_grid = 2*np.pi/(box_size)
            q_grid = np.sqrt(grid_coord[0]**2+grid_coord[1]**2+grid_coord[2]**2)*dq_grid
            
            dq = dq_grid
            nq = int(np.floor(dq_grid/dq*n_grid/2))
            qq = (np.arange(nq)+0.5)*dq
            index_q = np.floor(q_grid/dq) # q to the origin
            
            S_q = np.zeros(int(nq))
            
            for iq in range(int(nq)):
                #vq = 4*np.pi*(iq+0.5)**2/8
                #S_q[iq] = np.sum(S_q_lmn[index_q==iq])/vq/N
                S_q[iq] = np.average(S_q_lmn[index_q==iq])/N
                
        self.qq = qq
        self.S_q = S_q

#%% test
# # backbone
# # Coordinate of C atoms in each unit
# # unit_C = load('b_c.dat')';
# unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# # Degree of polymerization
# N_backbone = 10000

# # Chain stiffness
# a_backbone = 1e2

# # Unit persistence
# lambda_backbone = 1

# # call class
# chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
# tStart = time.time()
# chain01.chain()
# tEnd = time.time()
# print("It cost %f sec" % (tEnd - tStart))
# print('contour length = {:0.1f}'.format(chain01.l_contour))
# print('end-to-end distance = {:0.1f}'.format(chain01.l_end2end))
# chain01.plot()

#%%
# # call 'chain' function
# lc_backbone, Cc_backbone, O_backbone, n_backbone = chain_Rayleigh(N_backbone,a_backbone,lambda_backbone,unit_C)

# tEnd = time.time()
# print("It cost %f sec" % (tEnd - tStart))

# l_contour = np.sum(np.sqrt(np.sum(n_backbone**2,axis=0)))
# l_end2end = np.sqrt(np.sum((Cc_backbone[:,0]-Cc_backbone[:,-1])**2,axis=0))
# print('contour length = {:0.1f}'.format(l_contour))
# print('end-to-end length = {:0.1f}'.format(l_end2end))
       
# #%% plot
# plt.close('all')
# fig = plt.figure(figsize=(6, 6),Ni=192)
# ax = fig.add_subplot(projection='3d')

# ax.plot(Cc_backbone[0,:],Cc_backbone[1,:],Cc_backbone[2,:], 
#         '-', color='#D00000', linewidth=2, alpha = 0.75)
# # ax.plot(Cc_backbone[0,:],Cc_backbone[1,:],Cc_backbone[2,:], 
# #         'o', markeredgecolor='#800000', markerfacecolor='#D00000')

# ax.plot(Cc_backbone[0,0],Cc_backbone[1,0],Cc_backbone[2,0], 
#             'o', markeredgecolor='#800000', markerfacecolor='#D00000')
# ax.plot(Cc_backbone[0,-1],Cc_backbone[1,-1],Cc_backbone[2,-1], 
#             'o', markeredgecolor='#800000', markerfacecolor='#D00000')

# CM = np.mean(Cc_backbone,axis=1)
# CT = np.array([np.max(Cc_backbone[0,:])+np.min(Cc_backbone[0,:]),
#                np.max(Cc_backbone[1,:])+np.min(Cc_backbone[1,:]),
#                np.max(Cc_backbone[2,:])+np.min(Cc_backbone[2,:])])/2
# d_box = np.max([np.max(Cc_backbone[0,:])-np.min(Cc_backbone[0,:]),
#                 np.max(Cc_backbone[1,:])-np.min(Cc_backbone[1,:]),
#                 np.max(Cc_backbone[2,:])-np.min(Cc_backbone[2,:])])

# #ax.axis('off')
# ax.set_xlim([CT[0]-d_box/2, CT[0]+d_box/2])
# ax.set_ylim([CT[1]-d_box/2, CT[1]+d_box/2])
# ax.set_zlim([CT[2]-d_box/2, CT[2]+d_box/2])
# ax.set_box_aspect([1,1,1])
# # ax.set_xticklabels([])
# # ax.set_yticklabels([])
# # ax.set_zticklabels([])

# plt.show()

# #%%
# # fig2 = plt.figure(figsize=(6, 6),Ni=192)
# # ax2 = fig2.add_subplot()

# # ax2.plot(np.arange(N_backbone), np.sum(n_backbone**2,axis=0),'-')
# # ax2.set_yscale('log')
# # ax2.set_ylim([0.5, 2])

