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
filename = 'scatter_chain_block_0.mat'
scatter_dict = loadmat(filename)
qq_max = 1
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']

# load parameters
from scipy.io import loadmat
filename_stats = 'stats_block.mat'
filename_parameters = 'parameters_block.mat'
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
    filename = 'scatter_chain_block_{:d}.mat'.format(i)
    scatter_dict = loadmat(filename)
    S_q[:,i*np.shape(S_q_0)[1]:(i+1)*np.shape(S_q_0)[1]] = scatter_dict['S_q'][qq<qq_max,:]
    p[:,i*np.shape(p_0)[1]:(i+1)*np.shape(p_0)[1]] = scatter_dict['p']
    
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

qq = qq[qq<qq_max]
# S_q = S_q.T

F = S_q
F = (F[:,:].T*qq).T
F = np.log(F)
# F = F - np.mean(F,axis=0)
F = F.T
#%%prepare input
index_p_ra = (p[0] != set_ra[0])
index_p_a2 = (p[1] == set_a2[0])
index_p_f = (p[2] == set_f[0])
index_p = index_p_ra # bool
index_p = np.arange(len(p[1])) # all datapoints, bool

rng = np.random.default_rng(0)
index = np.arange(len(p[1]))[index_p]
rng.shuffle(index)
index_train = index[0:len(index)-100]
index_test = index[len(index)-100:]

# property to be presented 0
p_0 = p[0][index_train]
pc_0 = (p_0-np.nanmin(p_0))/(np.nanmax(p_0)-np.nanmin(p_0))

# property to be presented 1
p_1 = p[1][index_train]/1000
pc_1 = (p_1-np.nanmin(p_1))/(np.nanmax(p_1)-np.nanmin(p_1))

# property to be presented 2
p_2 = p[2][index_train]
pc_2 = (p_2-np.nanmin(p_2))/(np.nanmax(p_2)-np.nanmin(p_2))

# first moment
pm_1 = np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2

# second moment
pm_2 = (np.log(p_1)**2*(1-p_2) + np.log(p_0*p_1)**2*p_2 - 
    (np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2)**2)

# third moment
pm_3 = (p_2*(np.log(p_0*p_1)-pm_1)**3 + (1-p_2)*(np.log(p_1)-pm_1)**3)/np.sqrt(pm_2)**3

#%% GPR
X = F[index_train,:]
Y = pm_1

len_s = 0.254

sigma_y2 = 0.00416

kernel = RBF(len_s, (1e-3, 1e1)) + WhiteKernel(sigma_y2, (1e-4,1e-1))
gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, n_restarts_optimizer=10)
gp.fit(X, Y)

print("GPML kernel: %s" % gp.kernel_)
print("Log-marginal-likelihood: %.3f"
      % gp.log_marginal_likelihood(gp.kernel_.theta))

tFit = time.time()

#%% save model
import joblib
export_path_GPR = './saved_model/' 
model_name_GPR = 'model_GPR_mean'
export_name_GPR = export_path_GPR + model_name_GPR
joblib.dump(gp, export_name_GPR)
