# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 01:43:59 2022

@author: CHTUNG
"""

import numpy as np
import numpy.matlib
#from scipy.io import loadmat
#from scipy.io import savemat
import matplotlib.pyplot as plt
from WLM import WLChain
import time

# Degree of polymerization
N_backbone = 1000

# Chain stiffness
a_backbone = 1e8

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(N_backbone,a_backbone,lambda_backbone)
chain01.d_exc = 1

plt.close('all')
for i in range(1):
    tStart = time.time()
    chain01.d_exc = chain01.a*0.1*2
    chain01.apply_SA = 0
    chain01.chain()
    # chain01.chain_fix_val_free_rot()
    
    # # chain_grid method
    # chain01.d_exc = 1
    # chain01.kappa = 4
    # chain01.epsilon = 0.4
    # chain01.chain_grid()
    
    c_ave = chain01.cos_ave()
    d, corr = chain01.corr_o()
    
    # Z = chain01.Z
    # w,v = np.linalg.eig(Z)
    
    # iw = np.argsort(w)
    # w = w[iw]
    # v = v[:,iw]
    qq = (np.logspace(-4,0,65))
    chain01.scatter_direct(qq,n_merge=1)
    
    S_q = chain01.S_q
    
    tEnd = time.time()
    print("\'chain\' cost %f sec" % (tEnd - tStart))
    
    filename_ring = './figures/ring/ring_test_{:d}.png'.format(i+1)
    chain01.plot(filename=filename_ring, show_axes=0, save=0, end=1)

import Sk
S_q_rod = Sk.S_rod(qq,1000)

S_q_rod_bead = S_q

delta_S_q = S_q_rod/S_q

from scipy.io import savemat
filename = 'scatter_rod.mat'
mdic = {'S_q':S_q, 'qq':qq, 'delta_S_q':delta_S_q}
savemat(filename, mdic)

chain01.check_SA()