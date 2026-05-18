# -*- coding: utf-8 -*-
"""
Created on Wed May 20 15:05:20 2015

@author: jaekel
"""


########################
######## importing
########################
  
from __future__ import division

import numpy as np
from scipy import random

import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################

# number of observations lengths
N_max = 20
N = np.arange( 1, N_max+1 )

# number of trials per length
N_trial = int( 1e4 )

# parameter to be estimated and
# variance of additive noise
b = 1
sigma2 = .25

# initialize vector for storing variance of estimation error
mse_ML = np.zeros( N_max )
mse_MMSE = np.zeros( N_max )
mse_crb = sigma2 / N

# loop for lengths
for ind_n, val_n in enumerate(N):
    
    # initialize vector for estimation error of length n
    b_ML = np.zeros( N_trial )
    b_MMSE = np.zeros( N_trial )    
    
    # loop for trials
    for n_trial in np.arange(N_trial):
        y = b + random.normal( scale = np.sqrt(sigma2), size=(val_n,) )
        
        b_ML[ n_trial ] = np.average( y )
        b_MMSE[ n_trial ] = np.sum( y ) / ( val_n + sigma2 / b**2 )
        

    print( b, b_ML[ n_trial ], b_MMSE[ n_trial ] )

    mse_ML[ ind_n ] = np.var( b_ML )

    mse_MMSE[ ind_n ] = np.var( b_MMSE ) + ( b - np.mean(b_MMSE) )**2



# plotting
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

   
plt.figure()
plt.plot(mse_ML,'-o', label='$\widehat{\sigma^2}_{\mathrm{ML}}$')
plt.plot(mse_MMSE,'-o', label='$\hat{\sigma^2}_{\mathrm{MMSE}}$')
plt.plot(mse_crb, label='CRB')

plt.grid(True)
plt.legend(loc='upper right')
plt.xlabel('$N$')



plt.show()


############
print('Done!')