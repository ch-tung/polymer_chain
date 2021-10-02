# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 22:23:28 2021
Generate WLM chain trajectories
@author: CHTUNG
"""
import numpy as np
import numpy.matlib
from scipy.io import loadmat
from scipy.io import savemat

import time
tStart = time.time()

#%% define functions


#%% backbone
# Coordinate of C atoms in each unit
# unit_C = load('b_c.dat')';
unit_C = np.zeros((3,1)) # coordinate of C atoms in each unit

# Degree of polymerization
DP_backbone = 100

# Chain stiffness
a_backbone = 10

# Unit persistence
lambda_backbone = 1

DP = DP_backbone

n = np.zeros((3,DP))
l = np.zeros((3,DP))
lc = np.zeros((3,DP))
B = np.zeros((3,3))
C = np.zeros((3,3))
D = np.zeros((3,3))
R = np.zeros((3,3))
O = np.zeros((3,3,DP))

for i in range(DP):
    if i==0:
        n[:,i] = [1,0,0]
        l[:,i] = n[:,i]
        