# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 14:10:51 2021
Batch plot
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
import time
from WLM import WLChain
import matplotlib.pyplot as plt

n_plot = 20
plt.close('all')
for i in range(n_plot):
    #%% test
    # backbone
    # Coordinate of C atoms in each unit
    # unit_C = load('b_c.dat')';
    unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit
    
    # Degree of polymerization
    N_backbone = 5000
    
    # Chain stiffness
    a_backbone = np.array([2e2,2e1])
    
    # Unit persistence
    lambda_backbone = 1
    
    # call class
    chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
    tStart = time.time()
    chain01.d_exc = chain01.a*0.1*2
    chain01.apply_SA = 1
    
    # # chain_grid method
    # chain01.d_exc = 1
    # chain01.kappa = 5
    # chain01.epsilon = 0.1
    # chain01.chain_grid()
    
    chain01.f = 0.2
    chain01.chain_block()
    tEnd = time.time()
    print("It cost %f sec" % (tEnd - tStart))
    # print('contour length = {:0.1f}'.format(chain01.l_contour))
    # print('persistence length = {:0.1f}'.format(chain01.l_prstnc))
    # print('end-to-end distance = {:0.1f}'.format(chain01.l_end2end))
    # print('Rg = {:0.1f}'.format(chain01.Rg))
    filename_chain = './figures/chain/chain_block_{:d}.png'.format(i+1)
    chain01.plot_block(filename=filename_chain, show_axes=0, save=1, end=1)
    