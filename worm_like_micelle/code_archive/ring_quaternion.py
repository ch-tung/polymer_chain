# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 20:54:38 2021
Random walk polygon ring by quaternionic method
@author: CHTUNG
"""
import sys
sys.modules[__name__].__dict__.clear()

import numpy as np
import numpy.matlib
# import quaternion
#from scipy.io import loadmat
#from scipy.io import savemat
import matplotlib.pyplot as plt
from WLM import WLChain
import time

# parameters
chain01 = WLChain(N = 10000)

plt.close('all')
for i in range(1):
    tStart = time.time()
    #chain01.apply_SA = 0
    chain01.ring_q()
    tEnd = time.time()
    print("it cost %f sec" % (tEnd - tStart))
    
    filename_ring = './figures/ring/ring_test_{:d}.png'.format(i+1)
    chain01.plot(filename=filename_ring, show_axes=0, save=0, end=0)
