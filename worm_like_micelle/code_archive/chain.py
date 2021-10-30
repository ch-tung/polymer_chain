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

# Degree of polymerization
N_backbone = 10000

# Chain stiffness
a_backbone = 1e2

# Unit persistence
lambda_backbone = 1

# call class
chain01 = WLChain(N_backbone,a_backbone,lambda_backbone)
chain01.d_exc = 1

plt.close('all')
for i in range(1):
    tStart = time.time()
    #chain01.apply_SA = 0
    chain01.chain()
    tEnd = time.time()
    
    filename_ring = './figures/ring/ring_test_{:d}.png'.format(i+1)
    chain01.plot(filename=filename_ring, show_axes=0, save=0, end=0)

