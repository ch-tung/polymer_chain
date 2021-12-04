# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 22:15:38 2021
PCA on S(Q)
@author: CHTUNG
"""

import numpy as np
import numpy.matlib
from scipy.io import loadmat
# from scipy.io import savemat
# import matplotlib
import matplotlib.pyplot as plt
# from scipy import interpolate

#%% load
filename = 'scatter_chain_block.mat'
scatter_dict = loadmat(filename)
S_q  = scatter_dict['S_q']
p = scatter_dict['p']
set_a1 = set(p[0])
set_a2 = set(p[1])
set_f = set(p[2])

#%% SVD
F = S_q
# F = np.log(S_q)
F = F - np.mean(F,axis=0)
F = F
U, S, Vh = np.linalg.svd(F)

score_F = np.matmul(F.T,U)

#%% plot
index_p = (p[0] == list(set_a1)[0])
index_p = np.arange(64)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(projection='3d')
ax.scatter(score_F[index_p,0], score_F[index_p,1], score_F[index_p,2], 
           c=np.log(p[1][index_p]/p[0][index_p]),
           s=10)
ax.view_init(elev=25, azim=-135)
ax.set_xlabel('lv[0]')
ax.set_ylabel('lv[1]')
ax.set_zlabel('lv[2]')

plt.show()

