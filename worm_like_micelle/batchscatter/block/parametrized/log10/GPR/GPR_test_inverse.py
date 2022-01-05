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
filename = '../../../parametrized/scatter_chain_block_0.mat'
scatter_dict = loadmat(filename)
qq_max = 2
qq = scatter_dict['qq'][0,:]
S_q_0  = scatter_dict['S_q'][qq<qq_max,:]
p_0 = scatter_dict['p']

# load parameters
from scipy.io import loadmat
# filename_stats = 'stats_block.mat'
# filename_parameters = 'parameters_block.mat'
# stats_dict = loadmat(filename_stats)
# parameters_dict = loadmat(filename_parameters)

# parameters = parameters_dict['parameters']
# stats = stats_dict['statistics']

# set_stats0 = sorted(set(stats[:,0]))
# set_stats1 = sorted(set(stats[:,1]))
# set_stats2 = sorted(set(stats[:,2]))

n_mat =10
S_q = np.zeros((np.shape(S_q_0)[0],np.shape(S_q_0)[1]*n_mat))
p = np.zeros((np.shape(p_0)[0],np.shape(p_0)[1]*n_mat))
for i in range(n_mat):
    filename = '../../../parametrized/scatter_chain_block_{:d}.mat'.format(i)
    scatter_dict = loadmat(filename)
    S_q[:,i*np.shape(S_q_0)[1]:(i+1)*np.shape(S_q_0)[1]] = scatter_dict['S_q'][qq<qq_max,:]
    p[:,i*np.shape(p_0)[1]:(i+1)*np.shape(p_0)[1]] = scatter_dict['p']
    
set_ra = sorted(set(p[0]))
set_a2 = sorted(set(p[1]))
set_f = sorted(set(p[2]))

# fix the deviation high q limit
import Sk
q = qq
L = 1000
b = L*2
S_q_rod = Sk.S_rod(q,L)

d_th = 20
for i in range(S_q.shape[1]):
    chi = np.exp(-(q*d_th)**(-5))
    S_q_i = S_q[:,i]
    S_q[:,i] = chi*S_q_rod + (1-chi)*S_q_i
    # S_q_i[S_q_i<S_q_rod] = S_q_rod[S_q_i<S_q_rod]
    # S_q[:,i] = S_q_i

qq = qq[qq<qq_max]
# S_q = S_q.T

# SVD
def feature(S_q_in):
    F = S_q_in
    F = (F.T/S_q_rod.T).T
    # F = (F[:,:].T*qq).T
    F = np.log(F)
    # F = F - np.mean(F,axis=0)
    # F = np.gradient(F,axis=0)
    return F

F = feature(S_q)
F = F.T

#%%prepare input
index_p_ra = (p[0] != set_ra[0])
index_p_a2 = (p[1] == set_a2[0])
index_p_f = (p[2] == set_f[0])
index_p = index_p_ra # bool
# index_p = np.arange(len(p[1])) # all datapoints, bool

# rng = np.random.default_rng(0)
# index = np.arange(len(p[1]))[index_p]
# rng.shuffle(index)
# n_test = int(len(index)/10)
# index_train = index[0:len(index)-n_test]
# index_test = index[len(index)-n_test:]
index_test = index_p

# property to be presented 0
p_0 = p[0][index_test]
pc_0 = (p_0-np.nanmin(p_0))/(np.nanmax(p_0)-np.nanmin(p_0))

# property to be presented 1
p_1 = p[1][index_test]/1000
pc_1 = (p_1-np.nanmin(p_1))/(np.nanmax(p_1)-np.nanmin(p_1))

# property to be presented 2
p_2 = p[2][index_test]
pc_2 = (p_2-np.nanmin(p_2))/(np.nanmax(p_2)-np.nanmin(p_2))

# first moment
pm_1 = np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2

# second moment
pm_2 = (np.log(p_1)**2*(1-p_2) + np.log(p_0*p_1)**2*p_2 - 
    (np.log(p_1)*(1-p_2) + np.log(p_0*p_1)*p_2)**2)

# third moment
pm_3 = (p_2*(np.log(p_0*p_1)-pm_1)**3 + (1-p_2)*(np.log(p_1)-pm_1)**3)/np.sqrt(pm_2)**3

#%% GPR
X = F[index_test,:]
Y = [pm_1,pm_2,pm_3,np.log(p_1),np.log(p_0),p_2]

import joblib
export_path_GPR = './saved_model/' 
model_names = ['model_GPR_mean','model_GPR_var','model_GPR_skew','model_GPR_b1','model_GPR_rb','model_GPR_f']

plt.close('all')
gp_dict = {}
for i in range(6):
    export_name_GPR = export_path_GPR + model_names[i]
    gp_name_GPR = 'gp' + model_names[i]
    gp_dict[gp_name_GPR] = joblib.load(export_name_GPR)
    Y_predict = gp_dict[gp_name_GPR].predict(X, return_std=True)


    #%% plot prediction
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot()
    
    ax.errorbar(Y[i],Y_predict[0],yerr=Y_predict[1],fmt='.k',capsize=0)
    ax.scatter(Y[i],Y_predict[0],c=Y[i])
    ymin = min(Y[i])
    ymax = max(Y[i])
    ypmin = min(Y_predict[0])
    ypmax = max(Y_predict[0])
    ymin_all = min([ymin,ypmin])
    ymax_all = max([ymax,ypmax])
    dy_minmax = ymax_all-ymin_all
    ax.plot([ymin_all-dy_minmax/4,ymax_all+dy_minmax/4],[ymin_all-dy_minmax/4,ymax_all+dy_minmax/4],'-k')
    
    ax.set_xlabel('Y')
    ax.set_ylabel('Y_predict')