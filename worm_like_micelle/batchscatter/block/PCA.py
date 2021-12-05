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
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

#%% SVD
F = S_q
F = np.log(S_q)
F = F - np.mean(F,axis=0)
F = F
U, S, Vh = np.linalg.svd(F)

score_F = np.matmul(F.T,U)

#%% plot
index_p_a2 = (p[1] == set_a2[0])
index_p_ra = (p[0] == set_ra[2])
index_p_f = (p[2] == set_f[6])
index_p = index_p_a2
# index_p = np.arange(256)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(projection='3d')
ax.scatter(score_F[index_p,0], score_F[index_p,1], score_F[index_p,2], 
           c=np.log(p[0][index_p]),
           s=10)
ax.view_init(elev=25, azim=-135)
ax.set_xlabel('lv[0]')
ax.set_ylabel('lv[1]')
ax.set_zlabel('lv[2]')

plt.show()

