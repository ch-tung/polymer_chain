# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 14:53:21 2021
GPR on block polymers
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
filename = '../scatter_chain_block_0.mat'
scatter_dict = loadmat(filename)
qq_max = 1
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']

# load parameters
from scipy.io import loadmat
filename_stats = '../stats_block.mat'
filename_parameters = '../parameters_block.mat'
stats_dict = loadmat(filename_stats)
parameters_dict = loadmat(filename_parameters)

parameters = parameters_dict['parameters']
stats = stats_dict['statistics']

set_stats0 = sorted(set(stats[:,0]))
set_stats1 = sorted(set(stats[:,1]))
set_stats2 = sorted(set(stats[:,2]))

n_mat = 11
S_q = np.zeros((np.shape(S_q_0)[0],np.shape(S_q_0)[1]*n_mat))
p = np.zeros((np.shape(p_0)[0],np.shape(p_0)[1]*n_mat))
for i in range(n_mat):
    filename = '../scatter_chain_block_{:d}.mat'.format(i)
    scatter_dict = loadmat(filename)
    S_q[:,i*np.shape(S_q_0)[1]:(i+1)*np.shape(S_q_0)[1]] = scatter_dict['S_q'][qq<qq_max,:]
    p[:,i*np.shape(p_0)[1]:(i+1)*np.shape(p_0)[1]] = scatter_dict['p']
    
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

qq = qq[qq<qq_max]
S_q = S_q.T
#%%prepare input
index_p_ra = (p[0] != set_ra[0])
index_p_a2 = (p[1] == set_a2[0])
index_p_f = (p[2] == set_f[0])
index_p = index_p_ra
index_p = np.arange(len(p[1])) # all datapoints

# property to be presented 0
p_0 = p[0][index_p]
pc_0 = (p_0-np.nanmin(p_0))/(np.nanmax(p_0)-np.nanmin(p_0))

# property to be presented 1
p_1 = p[1][index_p]/1000
pc_1 = (p_1-np.nanmin(p_1))/(np.nanmax(p_1)-np.nanmin(p_1))

# property to be presented 2
p_2 = p[2][index_p]
pc_2 = (p_2-np.nanmin(p_2))/(np.nanmax(p_2)-np.nanmin(p_2))

# first moment
pm_1 = np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2

# second moment
pm_2 = (np.log(p_1)**2*(1-p_2) + np.log(p_0*p_1)**2*p_2 - 
    (np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2)**2)

# third moment
pm_3 = (p_2*(np.log(p_0*p_1)-pm_1)**3 + (1-p_2)*(np.log(p_1)-pm_1)**3)/np.sqrt(pm_2)**3

#%% GPR
X = S_q
Y = pm_1

kernel = RBF(1) + WhiteKernel(1)
gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, optimizer=None).fit(X, Y)

theta0 = np.logspace(-3, 0, 50) # len_s
theta1 = np.logspace(-4, -1, 50) # noise

print('start calculating LML')
tStart = time.time()

Theta0, Theta1 = np.meshgrid(theta0, theta1)
LML = [[gp.log_marginal_likelihood(np.log([Theta0[i, j], Theta1[i, j]]))
        for i in range(Theta0.shape[0])] for j in range(Theta0.shape[1])]
LML = np.array(LML).T

from scipy.io import savemat
mdic = {'len_s':theta0, 'noise':theta1, 'LML':LML}
savemat("LML_mean.mat", mdic)

tEnd = time.time()
print('job end')

print("It cost %f sec" % (tEnd - tStart))