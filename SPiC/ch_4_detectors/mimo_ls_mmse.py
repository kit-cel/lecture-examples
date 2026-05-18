# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:55:37 2015

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

import platform
import sys

if platform.system() == "Windows":
    sys.path.append('Z:\\Simulationen_Programme\kommunikationssysteme')
    
if platform.system() == "Linux":   
    sys.path.append('/home/jaekel/INT/Simulationen_Programme/kommunikationssysteme/')
    sys.path.append('/home/jaekel/Dokumente/Simulationen_Programme/kommunikationssysteme/')
        
import utils_communications as uc

########################
# main function
########################

# parameters of the MIMO system
N_t = 2
N_r = 2

# parameters for signal constellations and snr
M = 4
mod_scheme = 'QAM'
    
# find constellation points in the IQ plane
constellation = uc.find_constellation(M, mod_scheme)


# EbN0 to be used for the simulation
EbN0_db_min = 0
EbN0_db_max = 70    
EbN0_db_step = 10

EbN0_db_range = np.arange(EbN0_db_min, EbN0_db_max+EbN0_db_step, EbN0_db_step)

# map EbN0 to EsN0 for simulation    
EsN0_db_range = EbN0_db_range + 10*np.log10( np.log2(M) ) if M>2 else EbN0_db_range
EsN0_range =  10**(EsN0_db_range/10.)

# parameters limiting simulation time
max_errors = 1e2
max_syms = 1e6

# vectors of ser
ser_LS = np.zeros(len(EsN0_db_range))      
ser_MMSE = np.zeros(len(EsN0_db_range))      
      

# loop for snr        
for n_snr, snr_db in enumerate(EsN0_db_range):        
    
    # reset counters
    num_errors_LS = 0
    num_errors_MMSE = 0    
    num_syms = 0

    # loop for errors
    while ( num_syms < max_syms and num_errors_LS < max_errors):

        # generate data to be transmitted
        d = np.random.randint( 0, M, size=(N_t,) )
        s = uc.modulate(d, constellation) 
        
        # find channel matrix
        H = 1/np.sqrt(2) * ( np.random.randn(N_r, N_t) + 1j*np.random.randn(N_r, N_t) ) 
        x = np.dot( H, s)
        r = uc.add_awgn(x, snr_db, M, mod_scheme) 
        
        # demodule rx signal and look for errors
        H_LS = np.linalg.pinv(H)
        r_LS = np.dot(H_LS, r)
        
        H_MMSE = np.dot( np.linalg.pinv( np.dot( np.transpose(H.conj()), H) + 10**(-snr_db/10)*np.eye(N_t) ), np.transpose(H.conj()) )
        r_MMSE = np.dot( H_MMSE, r)
        
        d_est_LS = uc.demod_ML_awgn( r_LS, constellation) 
        d_est_MMSE = uc.demod_ML_awgn( r_MMSE, constellation) 
        
        # find errors
        for _k in np.arange(N_t):
            if abs(d_est_LS[_k]-d[_k])>10*np.spacing(1):
                num_errors_LS += 1
            if abs(d_est_MMSE[_k]-d[_k])>10*np.spacing(1):
                num_errors_MMSE += 1                
    
        num_syms += 1                    

    # find error rates
    ser_LS[n_snr] = num_errors_LS / num_syms
    ser_MMSE[n_snr] = num_errors_MMSE / num_syms    
    
    print( 'Es/N0 (dB) = %2.2f' %snr_db )


    


        
# theoretical curve
ser_theo = uc.ser_awgn_theory( np.array(EbN0_db_range), M, mod_scheme)
ser_fade_theo = 1./2. * ( 1- np.sqrt( EsN0_range/2 /(1+EsN0_range/2)))

# plotting
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

   
plt.figure()
plt.plot(EbN0_db_range, ser_fade_theo, label='Theo. '+mod_scheme)
plt.plot(EbN0_db_range, ser_LS, label='LS')
plt.plot(EbN0_db_range, ser_MMSE, label='MMSE')
plt.yscale('log')

plt.grid(True)
plt.legend(loc='lower left')
plt.xlabel('$N$')



plt.show()




############
print('Done!')
