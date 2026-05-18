# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 16:17:59 2014

@author: jaekel
"""

########################
# illustrating asymptotically unbiased estimation of psd
########################

# importing relevant stuff
import numpy as np
from scipy import signal
import own_utils_est_psd as own_spec

import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    # parameters: number of points in first resp. second function and in frequency domain
    N = 1e3
    M = N/5.0   
    
    # number of realizations for averaging    
    N_trials = 1e2

    # indices for acf
    N_acf = np.arange(-(N-1), N)
    
    # define windows
    rect = np.zeros(N)
    rect[ ((N+1)/2-M) : ((N+1)/2+M) ] = 1
    
    tria = 1.0-np.abs(N_acf)/M
    tria = map(lambda t: (t>0)*t, tria)
    
    hann = np.append( np.zeros(N-M), np.append(signal.hann(2*M), np.zeros( N-1-M)))
    hamming = np.append( np.zeros(N-M), np.append(signal.hamming(2*M), np.zeros( N-1-M)))    
    blackman = np.append( np.zeros( N-M), np.append(signal.blackman(2*M), np.zeros( N-1-M)))    
    
    
    # number of freq. points and freq. range
    N_freq = 512            
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
  
    # initialize two-dimensional arrays
    psd_c_tria = np.empty([N_trials, N_freq])
    psd_c_hann = np.empty([N_trials, N_freq])
    psd_c_hamming = np.empty([N_trials, N_freq])
    psd_c_blackman = np.empty([N_trials, N_freq])
            
    # loop for realizations
    n=0
    while n < N_trials:
    
        # wgn
        f_1 = np.sqrt(2)*np.random.normal(0.0, 1.0, N)
        
        # exponential
        #Omega_0 = 1.0
        #f_1 = np.exp(1j*Omega_0*np.arange(0, N))

        # activate to have filtered noise
        if 1:
            
            # filter parameters
            cutoff_freq = 1.0/4.0
            
            ripple_db = 60                      # ripples and transition width of the filter
            width = 1 / 5.0
            
            N_filter, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter

            taps = signal.firwin( N_filter, cutoff=cutoff_freq,  window=('kaiser', beta))
            f_1 = signal.lfilter(taps, 1.0, f_1)        
        
        # find acf estimations        
        #acf = est_acf(f_1, 'biased') 
        acf = own_spec.est_acf(f_1, 'biased')         
        
        
        # windows acfs and find correlograms
        psd_c_tria[n, :] = own_spec.find_correlogram( acf * tria, Ome)
        psd_c_hann[n, :] = own_spec.find_correlogram( acf * hann, Ome)
        psd_c_hamming[n, :] = own_spec.find_correlogram( acf * hamming, Ome)
        psd_c_blackman[n, :] = own_spec.find_correlogram( acf * blackman, Ome)

   
        n += 1
    
        # show progress
        done = float(n)/N_trials*100.0
        print 'Done: %3.2f percent' % done
     
     
    # find means and variations 
    psd_c_tria_mean = psd_c_tria.mean(axis=0)    
    psd_c_hann_mean = psd_c_hann.mean(axis=0)        
    psd_c_hamming_mean = psd_c_hamming.mean(axis=0)    
    psd_c_blackman_mean = psd_c_blackman.mean(axis=0)    
    
    psd_c_tria_deviation = psd_c_tria.var(axis=0)      
    psd_c_hann_deviation = psd_c_hann.var(axis=0)          
    psd_c_hamming_deviation = psd_c_hamming.var(axis=0)      
    psd_c_blackman_deviation = psd_c_blackman.var(axis=0)      
     
     
    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    
    # plot psd1
    plt.figure(2)    
    
    plt.subplot(221)    
    plt.plot(Ome, psd_c_tria_mean)      
    plt.plot(Ome, psd_c_tria_mean - psd_c_tria_deviation)          
    plt.plot(Ome, psd_c_tria_mean + psd_c_tria_deviation)  
      
    plt.title('triangular')    
    plt.grid(True); 
    plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.axis([-4, 4, 0, 3.5])
   
   
    plt.subplot(222)    
    plt.plot(Ome, psd_c_hann_mean)      
    plt.plot(Ome, psd_c_hann_mean - psd_c_hann_deviation)          
    plt.plot(Ome, psd_c_hann_mean + psd_c_hann_deviation)  
      
    plt.title('Hann')    
    plt.grid(True); 
    plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.axis([-4, 4, 0, 3.5])
  
  
    plt.subplot(223)        
    plt.plot(Ome, psd_c_hamming_mean)      
    plt.plot(Ome, psd_c_hamming_mean - psd_c_hamming_deviation)          
    plt.plot(Ome, psd_c_hamming_mean + psd_c_hamming_deviation)  
      
    plt.title('Hamming')    
    plt.grid(True); 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.axis([-4, 4, 0, 3.5])
    
    
    plt.subplot(224)        
    plt.plot(Ome, psd_c_blackman_mean)      
    plt.plot(Ome, psd_c_blackman_mean - psd_c_blackman_deviation)          
    plt.plot(Ome, psd_c_blackman_mean + psd_c_blackman_deviation)  
      
    plt.title('Blackman')    
    plt.grid(True); 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.axis([-4, 4, 0, 3.5])
    
    plt.show()
        
    
     
        





    
########################
# make it executable
########################
if __name__ == "__main__":
    main()