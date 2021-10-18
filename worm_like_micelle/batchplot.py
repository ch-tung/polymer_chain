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

n_plot = 1
plt.close('all')
for i in range(n_plot):
    #%% test
    # backbone
    # Coordinate of C atoms in each unit
    # unit_C = load('b_c.dat')';
    unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit
    
    # Degree of polymerization
    N_backbone = 10000
    
    # Chain stiffness
    a_backbone = 1e3
    
    # Unit persistence
    lambda_backbone = 1
    
    # call class
    chain01 = WLChain(N_backbone,a_backbone,lambda_backbone,unit_C)
    tStart = time.time()
    chain01.d_exc = 1
    chain01.apply_SA = 1
    chain01.chain()
    tEnd = time.time()
    print("It cost %f sec" % (tEnd - tStart))
    print('contour length = {:0.1f}'.format(chain01.l_contour))
    print('persistence length = {:0.1f}'.format(chain01.l_prstnc))
    print('end-to-end distance = {:0.1f}'.format(chain01.l_end2end))
    filename_chain = './figures/chain/chain_test_{:d}.png'.format(i+1)
    chain01.plot(filename=filename_chain, show_axes=0, save=0, end=0)
    