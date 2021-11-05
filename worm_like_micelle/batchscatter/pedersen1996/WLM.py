# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 22:23:28 2021
Generate WLM chain trajectories
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
#import quaternion
#from scipy.io import loadmat
#from scipy.io import savemat
import matplotlib.pyplot as plt
from scipy import interpolate
#import time

#%% define functions
def rotation(O,a):
    # quaternion
    phi_q = 2*(np.random.rand(1)-0.5)*np.pi
    delta = 0
    theta_q = np.sqrt(-np.log(1-np.random.rand(1)*(1-delta))*2/a)/2
    # ----------------------------------------
    # theta = 2*theta_q = sqrt(-ln(1-X)/a)
    # where X is a random variable in [0,1]
    #       a is the persistence length
    # ----------------------------------------
    
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

def rotation_dihedral(O,a):
    # quaternion
    phi_q = 2*(np.random.rand(1)-0.5)*np.pi
    theta_q = np.arctan(1/np.sqrt(a/2))/2
    # ----------------------------------------
    # theta = 2*theta_q
    # 2*a = [1+cos(theta)]/[1-cos(theta)]
    # where a is the persistence length
    # ----------------------------------------
    
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
    
    R = Rq
    
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
    i_diameter = int(np.ceil(5/3*d_exc/lambda_seg))
       
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
                    while SA == 0:
                        n_retry += 1
                        
                        if n_retry > 100:
                            abort = 1
                            #print('abort')
                            break
                            
                        d2_uv_min = np.min(np.sum((l[:,:i-i_diameter+1].T-l[:,i].T)**2,axis=1))
                        # d1_uv_min = np.min(np.max(np.abs(l[:,:i-1].T-l[:,i].T),axis=1))
                        # print(d1_uv_min)
                        
                        if d2_uv_min<d2_exc:
                        # if d1_uv_min<d_exc:
                            #print('retry ({:d})'.format(n_retry+1))
                            # n_retry+=1
                            R = rotation(O[:,:,i-1],a)
                
                            O[:,:,i] = R@O[:,:,i-1]
                            # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                            n[:,i] = O[:,1,i].reshape((3))
                            # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                            l[:,i] = l[:,i-1] + n[:,i]
                        else:
                            #if n_retry!=0:
                            #    print('retry (end)')
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

def chain_fix_val_free_rot(N, a, lambda_seg, unit_C, apply_SA=1, d_exc=1):
    d2_exc = d_exc**2
    i_diameter = int(np.ceil(5/3*d_exc/lambda_seg))
       
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
                    while SA == 0:
                        n_retry += 1
                        
                        if n_retry > 100:
                            abort = 1
                            #print('abort')
                            break
                            
                        d2_uv_min = np.min(np.sum((l[:,:i-i_diameter+1].T-l[:,i].T)**2,axis=1))
                        # d1_uv_min = np.min(np.max(np.abs(l[:,:i-1].T-l[:,i].T),axis=1))
                        # print(d1_uv_min)
                        
                        if d2_uv_min<d2_exc:
                        # if d1_uv_min<d_exc:
                            #print('retry ({:d})'.format(n_retry+1))
                            # n_retry+=1
                            R = rotation_dihedral(O[:,:,i-1],a)
                
                            O[:,:,i] = R@O[:,:,i-1]
                            # O[:,:,i] = O[:,:,i]/np.sqrt(np.sum(O[:,:,i]**2,axis=0))
                            n[:,i] = O[:,1,i].reshape((3))
                            # n[:,i] = n[:,i]/np.sqrt(np.sum(n[:,i]**2))
                            l[:,i] = l[:,i-1] + n[:,i]
                        else:
                            #if n_retry!=0:
                            #    print('retry (end)')
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

def ring_harmonic(N,n_harmonics,sigma):
    c_ring = np.zeros((3,N+1))
    # c_ring_deriv = np.zeros((3,N+1))
    
    theta = np.arange(N+1)/N*2*np.pi
    
    for i in range(3):
        phi_i = 2*np.pi*np.random.rand(1)
        
        weight = np.exp(-(np.arange(n_harmonics)+1)**2/sigma**2/2)
        weight = weight/np.sqrt(np.sum(weight**2))
        coeff_c_i = np.random.rand(n_harmonics)*weight
        coeff_s_i = np.random.rand(n_harmonics)*weight
        
        harmonics_c_i = np.cos(np.outer(theta,(np.arange(n_harmonics)+1)) + phi_i)*coeff_c_i
        harmonics_s_i = np.sin(np.outer(theta,(np.arange(n_harmonics)+1)) + phi_i)*coeff_s_i
        
        harmonics_i = harmonics_c_i + harmonics_s_i
        # harmonics_i_deriv = harmonics_s_i*(np.arange(n_harmonics)+1) + -harmonics_c_i*(np.arange(n_harmonics)+1)
        
        c_ring[i,:] = np.sum(harmonics_i,axis=1)
        # c_ring_deriv[i,:] = np.sum(harmonics_i_deriv,axis=1)
    
    arc_segment = np.sqrt(np.sum((c_ring[:,:-1]-c_ring[:,1:])**2,axis=0))
    
    # arc_segment = np.sqrt(np.sum(c_ring_deriv**2,axis=0))
    arc_sum = np.sum(arc_segment)
    
    arc_cum = np.zeros(N+1)
    for i in range(N+1):
        if i==0:
            arc_cum[i] = 0
            
        else:
            arc_cum[i] = arc_cum[i-1] + arc_segment[i-1]
            
    f_arc = interpolate.interp1d(arc_cum, theta, kind='quadratic', fill_value='extrapolate')
    
    arc_seq = np.arange(N+1)/(N+1)*arc_sum
    
    # print(arc_cum)
    # print(arc_seq)
    
    theta_interpolate = f_arc(arc_seq)
    
    f_ring = interpolate.interp1d(theta, c_ring, kind='quadratic', fill_value='extrapolate')
    Cc = f_ring(theta_interpolate)/arc_sum*(N+1)
    
    return Cc

#def ring_q(N,lmbda):
#    list_u = []
#    list_v = []
#
#    for i in range(N):
#        U = np.random.rand(2)
#        u_re = np.sqrt(-np.log(U[0]))*np.cos(2*np.pi*U[1])
#        u_im = np.sqrt(-np.log(U[0]))*np.sin(2*np.pi*U[1])
#        u = np.quaternion(u_re,u_im,0,0)
#        list_u.append(u)
#        V = np.random.rand(2)
#        v_re = np.sqrt(-np.log(V[0]))*np.cos(2*np.pi*V[1])
#        v_im = np.sqrt(-np.log(V[0]))*np.sin(2*np.pi*V[1])
#        v = np.quaternion(v_re,v_im,0,0)
#        list_v.append(v)
#    
#    uu = np.array(list_u)
#    vv = np.array(list_v)
#    
#    vv2 = vv - uu*np.sum(np.conjugate(uu)*vv)/np.sum(np.conjugate(uu)*uu)
#    uu = uu/np.sqrt(np.sum(np.conjugate(uu)*uu))
#    vv2 = vv2/np.sqrt(np.sum(np.conjugate(vv2)*vv2))
#    
#    h = uu + vv2*np.quaternion(0,0,1,0)
#    
#    h_Hopf = np.conjugate(h)*np.quaternion(0,1,0,0)*h
#    
#    n = quaternion.as_float_array(h_Hopf)[:,1:]*N
#    
#    l = np.zeros((N+1,3))
#    for i in range(N+1):
#        if i==0:
#            l[i,:] = 0
#            
#        else:
#            l[i,:] = l[i-1,:] + n[i-1,:]
#    
#    Cc = l.T*lmbda
#    
#    return Cc
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
        
    def chain(self):
        """
        Call the chain function acd calculate particle trajectory in WL-chain.
        """
        
        # call 'chain_Rayleigh' function
        self.lc, self.Cc, self.O, self.n = chain_Rayleigh(self.N,self.a,self.lmbda,self.unit_C,
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
        
        # call 'chain_Rayleigh' function
        self.lc, self.Cc, self.O, self.n = chain_fix_val_free_rot(self.N,self.a,self.lmbda,self.unit_C,
                                                          apply_SA=self.apply_SA,d_exc=self.d_exc)
        self.l_contour = np.sum(np.sqrt(np.sum(self.n**2,axis=0)))
        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
        self.l_prstnc = self.lmbda/(1-(1/np.tanh(self.a)-1/self.a))
        Cc_centered = self.Cc.T-np.mean(self.Cc.T,axis=0)
        self.Rg = np.sqrt(np.trace(Cc_centered.T@Cc_centered/self.N))
        #self.l_prstnc = np.dot(self.n[:,0].T,self.lc[:,-1])
        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
        
    def ring(self,n_harmonics,sigma):
        """
        Call the chain function and calculate particle trajectory in WL-chain.
        
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
        
#    def ring_q(self):
#        """
#        Call the chain function and calculate particle trajectory in WL-chain.
#        
#        Uehara, E., Tanaka, R., Inoue, M., Hirose, F., & Deguchi, T. (2014). 
#        Mean-square radius of gyration and hydrodynamic radius for topological 
#        polymers evaluated through the quaternionic algorithm. 
#        Reactive and Functional Polymers, 80, 48-56.
#        
#        numpy-quaternion package required https://github.com/moble/quaternion
#        """
#        
#        # call 'ring_harmonics' function
#        self.Cc = ring_q(self.N,self.lmbda)
#        self.l_end2end = np.sqrt(np.sum((self.Cc[:,0]-self.Cc[:,-1])**2,axis=0))
#        self.box = np.vstack((np.min(self.Cc, axis=1), np.max(self.Cc, axis=1)))
    
    def plot(self, filename=[], show_axes=1, save=0, end=1):
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
            n_q: int
                number of q points
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
        
    def scatter_direct(self, qq, n_merge=1):
        """
        Calculate scattering function.
        
        Args:
            n_q: int
                number of q points
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
            
        print('{:d} beads used to calculating S(q)'.format(Cc_merge.shape[1]))

        # two-point correlation
        n_list = N_merge
        # r_jk = self.Cc.T.reshape(n_list,1,3) - self.Cc.T.reshape(1,n_list,3)
        r_jk = Cc_merge.T.reshape(n_list,1,3) - Cc_merge.T.reshape(1,n_list,3)
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
        
        S_q = S_q/N_merge**2
            
        self.qq = qq
        self.S_q = S_q