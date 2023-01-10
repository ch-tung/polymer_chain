# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 11:53:50 2022

@author: CHTUNG
"""

import numpy as np
import numpy.matlib
from scipy.io import loadmat
# from scipy.io import savemat
# import matplotlib
import matplotlib.pyplot as plt
# from scipy import interpolate

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

import time
tStart = time.time()

#%% load
# load mat files
# grep shape
filename = '../scatter_chain_fvfr_SA.mat'
scatter_dict = loadmat(filename)
qq_max = 1
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]

qq = (np.logspace(-4,0,65))

import Sk

q = qq
L = 1000
b = L*2
S_q_th = Sk.Sk(q,L,b)
S_q_rod = Sk.S_rod(q,L)

plt.close('all')
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot()
ax.plot(qq,S_q_th)
ax.plot(qq,S_q_rod,'-k')
# ax.plot(qq[qq*b>3.1],S_q_th[qq*b>3.1])

ax.set_xscale('log')
ax.set_yscale('log')
