# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 22:23:28 2021
Generate WLM chain trajectories
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
# import quaternion
#from scipy.io import loadmat
#from scipy.io import savemat
import matplotlib.pyplot as plt
# from scipy import interpolate
#import time

#%% define functions
import f_rotation
rotation = f_rotation.rotation
rotation_dihedral = f_rotation.rotation_dihedral

import f_chain
chain_Rayleigh = f_chain.chain_Rayleigh
chain_Rayleigh_woSA = f_chain.chain_Rayleigh_woSA
chain_fix_val_free_rot = f_chain.chain_fix_val_free_rot
chain_fix_val_free_rot_woSA = f_chain.chain_fix_val_free_rot_woSA
chain_grid = f_chain.chain_grid
chain_grid_shear = f_chain.chain_grid_shear
chain_grid_woSA = f_chain.chain_grid_woSA
chain_grid_shear_woSA = f_chain.chain_grid_shear_woSA
chain_Rayleigh_block = f_chain.chain_Rayleigh_block
chain_Rayleigh_block_woSA = f_chain.chain_Rayleigh_block_woSA

import f_ring
ring_harmonic = f_ring.ring_harmonic
# ring_q = f_ring.ring_q

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
    l_prstnc = []
    Rg = []
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
        self.kappa = 1
        self.epsilon = 0
        self.f = 0
        self.grid = 'SC'
        
    def chain(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_Rayleigh' function
        if self.apply_SA == 0:
            self.lc, self.Cc, self.O, self.n = chain_Rayleigh_woSA(self.N,self.a,self.lmbda,self.unit_C,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc)
        else:
            self.lc, self.Cc, self.O, self.n = chain_Rayleigh(self.N,self.a,self.lmbda,self.unit_C,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc)
            
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.l_prstnc = self.lmbda/(1-(1/np.tanh(self.a)-1/self.a))
        Cc_centered = self.Cc.T-np.mean(self.Cc.T,axis=0)
        self.Rg = np.sqrt(np.trace(Cc_centered.T@Cc_centered/self.N))
        #self.l_prstnc = np.dot(self.n[:,0].T,self.lc[:,-1])
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def chain_block(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_Rayleigh' function
        if self.apply_SA == 0:
            self.lc, self.Cc, self.O, self.n, self.N1 = chain_Rayleigh_block_woSA(self.N,self.a,self.f,self.lmbda,self.unit_C,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc)
        else:
            self.lc, self.Cc, self.O, self.n, self.N1 = chain_Rayleigh_block(self.N,self.a,self.f,self.lmbda,self.unit_C,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc)
            
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.l_prstnc = self.lmbda/(1-(1/np.tanh(self.a)-1/self.a))
        Cc_centered = self.Cc.T-np.mean(self.Cc.T,axis=0)
        self.Rg = np.sqrt(np.trace(Cc_centered.T@Cc_centered/self.N))
        #self.l_prstnc = np.dot(self.n[:,0].T,self.lc[:,-1])
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def chain_fix_val_free_rot(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_fix_val_free_rot' function
        if self.apply_SA == 0:
            self.lc, self.Cc, self.O, self.n = chain_fix_val_free_rot_woSA(self.N,self.a,self.lmbda,self.unit_C,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc)
        else:
            self.lc, self.Cc, self.O, self.n = chain_fix_val_free_rot(self.N,self.a,self.lmbda,self.unit_C,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc)
            
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.l_prstnc = self.lmbda/(1-(1/np.tanh(self.a)-1/self.a))
        Cc_centered = self.Cc.T-np.mean(self.Cc.T,axis=0)
        self.Rg = np.sqrt(np.trace(Cc_centered.T@Cc_centered/self.N))
        #self.l_prstnc = np.dot(self.n[:,0].T,self.lc[:,-1])
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def chain_grid(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_Rayleigh' function
        if self.apply_SA == 0:
            self.lc, self.Cc, self.n, self.Z = chain_grid_woSA(self.N,self.kappa,self.epsilon,self.lmbda,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc,grid=self.grid)
        else:
            self.lc, self.Cc, self.n, self.Z = chain_grid(self.N,self.kappa,self.epsilon,self.lmbda,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc,grid=self.grid)
            
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.l_prstnc = 0.25/np.exp(-self.kappa)*self.lmbda
        Cc_centered = self.Cc.T-np.mean(self.Cc.T,axis=0)
        self.Rg = np.sqrt(np.trace(Cc_centered.T@Cc_centered/self.N))
        #self.l_prstnc = np.dot(self.n[:,0].T,self.lc[:,-1])
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def chain_grid_shear(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_Rayleigh' function
        if self.apply_SA == 0:
            self.lc, self.Cc, self.n, self.Z = chain_grid_shear_woSA(self.N,self.kappa,self.epsilon,self.lmbda,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc,grid=self.grid)
        else:
            self.lc, self.Cc, self.n, self.Z = chain_grid_shear(self.N,self.kappa,self.epsilon,self.lmbda,
                                                              apply_SA=self.apply_SA,d_exc=self.d_exc,grid=self.grid)
            
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.l_prstnc = 0.25/np.exp(-self.kappa)*self.lmbda
        Cc_centered = self.Cc.T-np.mean(self.Cc.T,axis=0)
        self.Rg = np.sqrt(np.trace(Cc_centered.T@Cc_centered/self.N))
        #self.l_prstnc = np.dot(self.n[:,0].T,self.lc[:,-1])
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def cos_ave(self):
        c_ave = 0
        for i in range(self.N-1):
            c_ave = c_ave + np.dot(self.n[:,i],self.n[:,i+1])
            
        c_ave = c_ave/self.N 
        return c_ave
    
    def corr_o(self):
        d = np.ceil(10**(np.arange(64)/8)).astype('int')
        d = d[d<self.N]
        print(len(d))
        corr = np.zeros(len(d))
        for i, d_i in enumerate(d):
            print(i)
            for j in range(self.N-d_i):
                corr[i] = corr[i] + np.dot(self.n[:,j],self.n[:,j+d_i])
            corr[i] = corr[i]/(self.N-d_i)
            
        return d, corr
        
    def ring(self,n_harmonics,sigma):
        """
        Call the ring function and calculate particle trajectory in ring polymer.
        
        Args:
            n_harmonics: int
                harmonics used in fourier series
                
            sigma: float
                controlled the spread of k-distribution (p(k) = exp(-k^2/(2*sigma^2)))
        """
        
        # call 'ring_harmonics' function
        self.Cc = ring_harmonic(self.N,n_harmonics,sigma)
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    # def ring_q(self):
    #     """
    #     Call the ring function and calculate particle trajectory in ring polymer.
        
    #     Uehara, E., Tanaka, R., Inoue, M., Hirose, F., & Deguchi, T. (2014). 
    #     Mean-square radius of gyration and hydrodynamic radius for topological 
    #     polymers evaluated through the quaternionic algorithm. 
    #     Reactive and Functional Polymers, 80, 48-56.
        
    #     numpy-quaternion package required https://github.com/moble/quaternion
    #     """
        
    #     # call 'ring_q' function
    #     self.Cc = ring_q(self.N,self.lmbda)
    #     self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
    #     self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
    
    def plot(self, filename=[], show_axes=1, save=0, end=1, axeslabel='off'):
        """
        Plot polymer chain.
        
        Args:
            filename: str
                path of the generated figure
                
            show_axes: boolean
            
            save: boolean
            
            end: boolean
                whether to display the end-point of loop
        """
        
        #plt.close('all')
        fig = plt.figure(figsize=(6, 6),dpi=192)
        ax = fig.add_subplot(projection='3d')
        
        ax.plot(self.Cc[0,:],self.Cc[1,:],self.Cc[2,:], 
                '-', color='#D00000', linewidth=2, alpha = 0.75)
        # ax.plot(self.Cc[0,:],self.Cc[1,:],self.Cc[2,:], 
        #         'o', markeredgecolor='#800000', markerfacecolor='#D00000')
        
        # plot chain end
        if end==1:
            ax.plot(self.Cc[0,0],self.Cc[1,0],self.Cc[2,0], 
                        'o', markeredgecolor='#000080', markerfacecolor='#0000D0')
            ax.plot(self.Cc[0,-1],self.Cc[1,-1],self.Cc[2,-1], 
                        'o', markeredgecolor='#008000', markerfacecolor='#00D000')
        
        #CM = np.mean(Cc_backbone,axis=1)
        CT = np.array([np.max(self.Cc[0,:])+np.min(self.Cc[0,:]),
                       np.max(self.Cc[1,:])+np.min(self.Cc[1,:]),
                       np.max(self.Cc[2,:])+np.min(self.Cc[2,:])])/2
        d_box = np.max([np.max(self.Cc[0,:])-np.min(self.Cc[0,:]),
                        np.max(self.Cc[1,:])-np.min(self.Cc[1,:]),
                        np.max(self.Cc[2,:])-np.min(self.Cc[2,:])])
        
        if axeslabel=='on':
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
        
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
        
    def plot_block(self, filename=[], show_axes=1, save=0, end=1):
        """
        Plot polymer chain.
        
        Args:
            filename: str
                path of the generated figure
                
            show_axes: boolean
            
            save: boolean
            
            end: boolean
                whether to display the end-point of loop
        """
        
        #plt.close('all')
        fig = plt.figure(figsize=(6, 6),dpi=192)
        ax = fig.add_subplot(projection='3d')
        
        ax.plot(self.Cc[0,:self.N1],self.Cc[1,:self.N1],self.Cc[2,:self.N1], 
                '-', color='#D00000', linewidth=2, alpha = 0.75)
        ax.plot(self.Cc[0,self.N1:],self.Cc[1,self.N1:],self.Cc[2,self.N1:], 
                '-', color='#0000D0', linewidth=2, alpha = 0.75)
        # ax.plot(self.Cc[0,:],self.Cc[1,:],self.Cc[2,:], 
        #         'o', markeredgecolor='#800000', markerfacecolor='#D00000')
        
        # plot chain end
        if end==1:
            ax.plot(self.Cc[0,0],self.Cc[1,0],self.Cc[2,0], 
                        'o', markeredgecolor='#008000', markerfacecolor='#00D000')
            ax.plot(self.Cc[0,-1],self.Cc[1,-1],self.Cc[2,-1], 
                        'o', markeredgecolor='#008000', markerfacecolor='#00D000')
        
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
    
    def close(self):
        plt.close('all')
        
    def scatter_grid(self, n_grid=256, approx_1D=0, box_size=1e4):
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
        #box_size = N
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

    def scatter_grid_direct(self, qq, n_grid=128, box_size=1e4):
        """
        Calculate scattering function.
        
        Args:
            n_grid: int
                number of grid points
            qq: array
                wave vectors
            box_size: float
                size of cubic box containing the polymer chain
        """
        
        N = self.N
        chain_box = self.box
        
        # box_size = np.max(chain_box[1,:]-chain_box[0,:])+1
        # box_size = N
        grid_size = (box_size)/n_grid
        Cc_relative = self.Cc.T-chain_box[0,:] # relative position of WL-chain in the box
        bead_coord = np.floor(Cc_relative/grid_size).astype('int')
        
        # density in real space
        rho_r = np.zeros((n_grid,n_grid,n_grid))
        
        for i in range(N):
            rho_r[bead_coord[i,0],bead_coord[i,1],bead_coord[i,2]] += 1
        
        list_rho_r_all = rho_r.reshape(-1)
        index_rho_r = np.nonzero(list_rho_r_all)[0]
        list_rho_r = list_rho_r_all[index_rho_r]
        coord_rho_r_x = np.floor_divide(index_rho_r, n_grid**2)
        coord_rho_r_y = np.floor_divide(index_rho_r-coord_rho_r_x*n_grid**2, n_grid)
        coord_rho_r_z = index_rho_r-coord_rho_r_x*n_grid**2-coord_rho_r_y*n_grid
        coord_rho_r = np.vstack((coord_rho_r_x,coord_rho_r_y,coord_rho_r_z))*box_size/n_grid
        
        # two-point correlation
        n_list = len(list_rho_r)
        r_jk = coord_rho_r.T.reshape(n_list,1,3) - coord_rho_r.T.reshape(1,n_list,3)
        d_jk = np.sqrt(np.sum(r_jk**2,axis=2))
        n_jk = np.outer(list_rho_r, list_rho_r)
        rho_jk = n_jk/np.sum(n_jk)
        
        # radial average
        # dq_grid = 2*np.pi/(box_size)
        # dq = dq_grid
        # nq = int(np.floor(dq_grid/dq*n_grid/2))
        # qq0 = np.arange(nq)+0.5
        # qq = qq0*dq
        # qq0 = 2*np.pi/(np.logspace(1,5,n_q))
        nq = len(qq)
        # qq = qq0 
        
        S_q = np.zeros(int(nq))
        d_jk_list = d_jk[d_jk!=0]
        rho_jk_list = rho_jk[d_jk!=0]
        
        for iq in range(int(nq)):
            sinqr_qr = rho_jk_list*np.sin(qq[iq]*d_jk_list)/(qq[iq]*d_jk_list)
            S_q[iq] = np.sum(sinqr_qr[np.isnan(sinqr_qr)==0])
                
        self.qq = qq
        self.S_q = S_q
        
    def scatter_direct(self, qq, n_merge=1, p_sub=1.0):
        """
        Calculate scattering function.
        
        Args:
            qq: array
                wave vectors
            n_merge: int
                merge consecutive n_merge beads into one bead
        """
        
        N = self.N
        # chain_box = self.box
        
        # merge beads
        N_merge = int(N/n_merge)
        Cc_merge = np.zeros((3,N_merge))
        for i in range(N_merge):
            Cc_merge[:,i] = np.mean(self.Cc[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
            
#         print('{:d} beads used to calculate S(q)'.format(Cc_merge.shape[1]))

        # two-point correlation
        n_list = int(N_merge*p_sub)
        i_list = np.random.choice(np.arange(N_merge), size=n_list)
        # r_jk = self.Cc.T.reshape(n_list,1,3) - self.Cc.T.reshape(1,n_list,3)
        # r_jk = Cc_merge.T.reshape(n_list,1,3) - Cc_merge.T.reshape(1,n_list,3)
        r_jk = Cc_merge[:,i_list].T.reshape(n_list,1,3) - Cc_merge[:,i_list].T.reshape(1,n_list,3)
        d_jk = np.sqrt(np.sum(r_jk**2,axis=2))
        
        # radial average
        # dq_grid = 2*np.pi/(box_size)
        # dq = dq_grid
        # nq = int(np.floor(dq_grid/dq*n_grid/2))
        # qq0 = np.arange(nq)+0.5
        # qq = qq0*dq
        # qq0 = 2*np.pi/(np.logspace(1,5,n_q))
        nq = len(qq)
        # qq = qq0 
        
        S_q = np.zeros(int(nq))
        d_jk_list = d_jk[d_jk!=0]
        
        for iq in range(int(nq)):
            sinqr_qr = np.sin(qq[iq]*d_jk_list)/(qq[iq]*d_jk_list)
            S_q[iq] = np.sum(sinqr_qr[np.isnan(sinqr_qr)==0])
        
        S_q = S_q/n_list**2
            
        self.qq = qq
        self.S_q = S_q
        
    def scatter_direct_pw(self, qq, n_merge=1):
        """
        Calculate scattering function by taking average at the stage of 
        plane waves superposition.
        *Applys only for isotropic systems*
        
        Args:
            qq: array
                wave vectors
            n_merge: int
                merge consecutive n_merge beads into one bead
        """
        
        N = self.N
        # chain_box = self.box
        
        # merge beads
        N_merge = int(N/n_merge)
        Cc_merge = np.zeros((3,N_merge))
        for i in range(N_merge):
            Cc_merge[:,i] = np.mean(self.Cc[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
            
        print('{:d} beads used to calculate S(q)'.format(Cc_merge.shape[1]))
        
        phi = np.zeros((3,len(qq)),dtype = 'complex_')
        for i in range(len(qq)):
            phi[:,i] = np.sum(np.exp((1.j)*qq[i]*Cc_merge),axis=1)
        
        phi = phi
        
        def abs2(x):
            return x.real**2 + x.imag**2
        
        S_q = np.sum(abs2(phi),axis=0)/3/N_merge**2
            
        self.qq = qq
        self.S_q = S_q
        
    def scatter_direct_block(self, qq, n_merge=1):
        """
        Calculate scattering function.
        
        Args:
            qq: array
                wave vectors
            n_merge: int
                merge consecutive n_merge beads into one bead
        """
        
        N = self.N
        f = self.f
        N1 = int(N*f)
        N2 = N-N1
        # chain_box = self.box
        
        # merge beads
        N_merge = int(N/n_merge)
        # Cc_merge = np.zeros((3,N_merge))
        # for i in range(N_merge):
        #     Cc_merge[:,i] = np.mean(self.Cc[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
        
        Cc_1 = self.Cc[:,:self.N1]
        Cc_2 = self.Cc[:,self.N1:]
        
        # # merge beads 12
        N1_merge = int(N1/n_merge)
        N2_merge = int(N2/n_merge)
        Cc_1_merge = np.zeros((3,N1_merge))
        Cc_2_merge = np.zeros((3,N2_merge))
        for i in range(N1_merge):
            Cc_1_merge[:,i] = np.mean(Cc_1[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
        for i in range(N2_merge):
            Cc_2_merge[:,i] = np.mean(Cc_2[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
            
        print('{:d} + {:d} beads used to calculate S(q)'.format(N1_merge,N2_merge))

        # two-point correlation              
        nq = len(qq)
        
        #S_q_11
        r_11_jk = Cc_1_merge.T.reshape(N1_merge,1,3) - Cc_1_merge.T.reshape(1,N1_merge,3)
        d_11_jk = np.sqrt(np.sum(r_11_jk**2,axis=2))
        
        S_q_11 = np.zeros(int(nq))
        d_11_jk_list = d_11_jk[d_11_jk!=0]

        #S_q_22
        r_22_jk = Cc_2_merge.T.reshape(N2_merge,1,3) - Cc_2_merge.T.reshape(1,N2_merge,3)
        d_22_jk = np.sqrt(np.sum(r_22_jk**2,axis=2))
        
        S_q_22 = np.zeros(int(nq))
        d_22_jk_list = d_22_jk[d_22_jk!=0]

        #S_q_12
        r_12_jk = Cc_1_merge.T.reshape(N1_merge,1,3) - Cc_2_merge.T.reshape(1,N2_merge,3)
        d_12_jk = np.sqrt(np.sum(r_12_jk**2,axis=2))
        
        S_q_12 = np.zeros(int(nq))
        d_12_jk_list = d_12_jk[d_12_jk!=0]
        
        for iq in range(int(nq)):
            sinqr_qr_11 = np.sin(qq[iq]*d_11_jk_list)/(qq[iq]*d_11_jk_list)
            S_q_11[iq] = np.sum(sinqr_qr_11[np.isnan(sinqr_qr_11)==0])
            sinqr_qr_12 = np.sin(qq[iq]*d_12_jk_list)/(qq[iq]*d_12_jk_list)
            S_q_12[iq] = np.sum(sinqr_qr_12[np.isnan(sinqr_qr_12)==0])
            sinqr_qr_22 = np.sin(qq[iq]*d_22_jk_list)/(qq[iq]*d_22_jk_list)
            S_q_22[iq] = np.sum(sinqr_qr_22[np.isnan(sinqr_qr_22)==0])
                
        S_q_11 = S_q_11/N_merge**2        
        S_q_22 = S_q_22/N_merge**2
        S_q_12 = S_q_12/N_merge**2*2
            
        self.qq = qq
        self.S_q_11 = S_q_11
        self.S_q_22 = S_q_22
        self.S_q_12 = S_q_12
        self.S_q = S_q_12+S_q_11+S_q_22
        
    def scatter_direct_aniso(self, qq, n_merge=1):
        """
        Calculate 2D spectra.
        
        Args:
            qq: array
                wave vectors
            n_merge: int
                merge consecutive n_merge beads into one bead
        """
        
        N = self.N
        # chain_box = self.box
        
        # merge beads
        N_merge = int(N/n_merge)
        Cc_merge = np.zeros((3,N_merge))
        for i in range(N_merge):
            Cc_merge[:,i] = np.mean(self.Cc[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
            
        print('{:d} beads used to calculate S(q)'.format(Cc_merge.shape[1]))

        # two-point correlation
        n_list = N_merge
        # r_jk = self.Cc.T.reshape(n_list,1,3) - self.Cc.T.reshape(1,n_list,3)
        r_jk = Cc_merge.T.reshape(n_list,1,3) - Cc_merge.T.reshape(1,n_list,3)
        d_jk = np.sqrt(np.sum(r_jk**2,axis=2))
        d_jk_list = d_jk[d_jk!=0]
        r_jk_list = np.zeros((len(d_jk_list),3))
        for i in range(3):
            r_jk_i = r_jk[:,:,i]
            r_jk_list[:,i] = r_jk_i[d_jk!=0]
        
        nq = len(qq)
        
        '''
        2D spectrum along 
        velocity gradient–vorticity(y–z), flow–vorticity (x–z), and flow–velocity gradient (x–y) planes
        '''
        
        def abs2(x):
            return x.real**2 + x.imag**2
        
        S_q_2D = np.zeros((int(nq)*2+1,int(nq)*2+1,3))
        qq_2D = np.concatenate((-np.flip(qq), np.array([0.0]), qq))
        i_axes_list = [[1,2],[0,2],[0,1]]
        for i, i_axes in enumerate(i_axes_list):
            print('{:01d}{:01d} plane'.format(i_axes[0]+1,i_axes[1]+1))
            for iqx in range(len(qq_2D)):
                qqx = qq_2D[iqx]
                for iqy in range(len(qq)+1):
                    qqy = qq_2D[iqy]
                    # if (qqx*2+qqy**2) != 0:
                    qr_xy = qqx*r_jk_list[:,i_axes[0]] + qqy*r_jk_list[:,i_axes[1]]
                    # sinqr_qr_2D = np.sin(qr_xy)/(qr_xy)
                    # S_q_2D[iqx,iqy,i] = np.sum(sinqr_qr_2D[np.isnan(sinqr_qr_2D)==0])
                    phi = np.exp((-1.j)*qr_xy)
                    S_q_2D[iqx,iqy,i] = abs2(np.mean(phi))
            
            # S(q) = S(-q)
            for iqx in range(len(qq_2D)):
                for iqy in range(len(qq_2D)):
                    qqy = qq_2D[iqy]
                    if qqy>0:
                        S_q_2D[iqx,iqy,i] = S_q_2D[len(qq_2D)-1-iqx,len(qq_2D)-1-iqy,i]
                        
            # S_q_2D[:,:,i] = S_q_2D[:,:,i]/N_merge**2
            # S_q_2D[len(qq),len(qq),i] = 1
            
        self.qq_2D = qq_2D
        self.S_q_2D = S_q_2D
            
                                    
        '''
        1D spectrum
        '''
        S_q = np.zeros(int(nq))
        for iq in range(int(nq)):
            sinqr_qr = np.sin(qq[iq]*d_jk_list)/(qq[iq]*d_jk_list)
            S_q[iq] = np.sum(sinqr_qr[np.isnan(sinqr_qr)==0])
        
        S_q = S_q/N_merge**2
            
        self.qq = qq
        self.S_q = S_q
        
    def scatter_direct_RSHE(self, qq, rr=[], n_merge=1, calculate_g_r=0):
        """
        Calculate scattering function.
        
        Args:
            qq: array
                wave vectors
            rr: array
                pair distances
            n_merge: int
                merge consecutive n_merge beads into one bead
            calculate_g_r: 0 or 1
                if 1, calculate the RSHE of real space correlations
        """
        
        N = self.N
        # chain_box = self.box
        
        # merge beads
        N_merge = int(N/n_merge)
        Cc_merge = np.zeros((3,N_merge))
        for i in range(N_merge):
            Cc_merge[:,i] = np.mean(self.Cc[:,i*n_merge:(i*n_merge+n_merge)],axis=1)
            
        print('{:d} beads used to calculate S(q)'.format(Cc_merge.shape[1]))

        # two-point correlation
        n_list = N_merge
        # r_jk = self.Cc.T.reshape(n_list,1,3) - self.Cc.T.reshape(1,n_list,3)
        r_jk = Cc_merge.T.reshape(n_list,1,3) - Cc_merge.T.reshape(1,n_list,3)
        d_jk = np.sqrt(np.sum(r_jk**2,axis=2))
        d_jk_list = d_jk[d_jk!=0]
        r_jk_list = np.zeros((len(d_jk_list),3))
        for i in range(3):
            r_jk_i = r_jk[:,:,i]
            r_jk_list[:,i] = r_jk_i[d_jk!=0]
        
        nq = len(qq)
        nr = len(rr)
        
        S_q = np.zeros(int(nq))
        S_q_lm = np.zeros((int(nq),6))
        g_r = np.zeros(int(nr))
        g_r_lm = np.zeros((int(nr),6))
        
        RSHE_coeff = [np.sqrt(1/np.pi)/2,
                      np.sqrt(15/np.pi)/2,np.sqrt(15/np.pi)/2,np.sqrt(5/np.pi)/4,np.sqrt(15/np.pi)/2,np.sqrt(15/np.pi)/4]/(np.sqrt(1/np.pi)/2)
        
        def j0(x):
            return np.sin(x)/(x)
        
        def j2(x):
            return np.sin(x)/(x)*(3/x**2-1)-3*np.cos(x)/x**2        
        
        for iq in range(int(nq)):
            j0_qr = j0(qq[iq]*d_jk_list)
            j2_qr = j2(qq[iq]*d_jk_list)
            S_q[iq] = np.sum(j0_qr[np.isnan(j0_qr)==0])
            
            Y0mq = RSHE_coeff[0]
            S_q_lm[iq,0] = np.sum((j0_qr*Y0mq)[np.isnan(j0_qr)==0])
            
            Y2mq = [RSHE_coeff[1]*r_jk_list[:,0]*r_jk_list[:,1]/d_jk_list**2,
                    RSHE_coeff[2]*r_jk_list[:,1]*r_jk_list[:,2]/d_jk_list**2,
                    RSHE_coeff[3]*(2*r_jk_list[:,2]**2-r_jk_list[:,0]**2-r_jk_list[:,1]**2)/d_jk_list**2,
                    RSHE_coeff[4]*r_jk_list[:,0]*r_jk_list[:,2]/d_jk_list**2,
                    RSHE_coeff[5]*(r_jk_list[:,0]**2-r_jk_list[:,1]**2)/d_jk_list**2]
            for im in range(5):
                S_q_lm[iq,im+1] = np.sum((j2_qr*Y2mq[im])[np.isnan(j2_qr)==0])
        
        S_q = S_q/N_merge**2
        S_q_lm = S_q_lm/N_merge**2
        
        self.qq = qq
        self.S_q = S_q
        self.S_q_lm = S_q_lm
        
        if calculate_g_r == 1:
            for ir in range(int(nr)):
                if ir == 0:
                    continue
                index_r = (d_jk_list>rr[ir-1])&(d_jk_list<rr[ir])
                n_r = np.sum(index_r)
                d_r = rr[ir]-rr[ir-1]
                g_r[ir] = n_r/d_r
                
                g_r_lm[ir,0] = n_r/d_r
                
                Y2mq = [RSHE_coeff[1]*r_jk_list[index_r,0]*r_jk_list[index_r,1]/d_jk_list[index_r]**2,
                        RSHE_coeff[2]*r_jk_list[index_r,1]*r_jk_list[index_r,2]/d_jk_list[index_r]**2,
                        RSHE_coeff[3]*(2*r_jk_list[index_r,2]**2-r_jk_list[index_r,0]**2-r_jk_list[index_r,1]**2)/d_jk_list[index_r]**2,
                        RSHE_coeff[4]*r_jk_list[index_r,0]*r_jk_list[index_r,2]/d_jk_list[index_r]**2,
                        RSHE_coeff[5]*(r_jk_list[index_r,0]**2-r_jk_list[index_r,1]**2)/d_jk_list[index_r]**2]
                for im in range(5):
                    g_r_lm[ir,im+1] = np.sum(Y2mq[im])/d_r
            
            self.rr = rr
            self.g_r = g_r
            self.g_r_lm = g_r_lm
            
    def segment_corr(self, rr=[]):
        """
        Calculate scattering function.
        
        Args:
            rr: array
                pair distances
        """
        
        N = self.N
        segment = self.n # segment vector

        # two-point correlation
        n_list = N
        r_jk = self.Cc.T.reshape(n_list,1,3) - self.Cc.T.reshape(1,n_list,3)
        d_jk = np.sqrt(np.sum(r_jk**2,axis=2))

        nr = len(rr)
        
        M_jk = np.zeros((3,3,n_list,n_list))
        for l in range(3):
            for m in range(3):
                M_jk[l,m,:,:] = np.outer(segment[l,:],segment[m,:])
        M_jk_lm = np.reshape(M_jk,(3,3,n_list**2))
        
        cos_lm = np.trace(M_jk_lm, axis1=0, axis2=1)
        P2_lm = (3*cos_lm**2-1)/2

        cos_r = np.zeros(int(nr))
        P2_r = np.zeros(int(nr))
        M_corr = np.zeros((3,3,nr))
        
        dr = rr[1]-rr[0]
        index_jk = np.floor(d_jk/dr)

        for ir in range(nr):
            index_r = np.reshape(index_jk==ir,n_list**2)
            cos_r[ir] = np.mean(cos_lm[index_r])
            P2_r[ir] = np.mean(P2_lm[index_r])
            M_corr[:,:,ir] = np.mean(M_jk_lm[:,:,index_r],axis=2)
        
        self.rr = rr
        self.cos_r = cos_r
        self.M_corr = M_corr
        self.P2_r = P2_r
        
    def check_SA(self):
        """
        check self avoiding
        """
        d2_exc = self.d_exc**2
        i_diameter = (np.ceil(5/3*self.d_exc/self.lmbda)).astype(int)
        n_intersection = 0
        for i in range(i_diameter,self.N):
            d2_ij = np.min(np.sum((self.Cc[:,i].T-self.Cc[:,:i-i_diameter+1].T)**2,axis=1))
            if d2_ij < d2_exc:
                n_intersection += 1
                # break
                    
        if n_intersection != 0:
            print('Self intersection detected! ({:d})'.format(n_intersection))
        else:
            print('No self intersection detected')
        