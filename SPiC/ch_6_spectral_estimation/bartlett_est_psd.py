# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 08:20:46 2014

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
    N_vec = np.arange(0, N)

    # member per group for Bartlett's method
    M = N/10.0
       
    # number of freq. points and freq. range
    N_freq = 512            
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
    # number of realizations for averaging    
    N_trials = 1e2
    
    # initializa arrays    
    psd_p = np.empty([N_trials, N_freq])
    psd_bartlett = np.empty([N_trials, N_freq])


    # loop for realizations, averaging and variance
    n=0
    while n < N_trials:
    
        # generate signal
        f_1 = np.sqrt(2) * np.random.normal(0.0, 1.0, N)
        
        Omega_0 = 1.0
        Omega_1 = 1.2
        f_1 = np.sin(Omega_0 * N_vec) + np.sin(Omega_1 * N_vec)  + np.random.normal(0.0, 1.0, size = N)

        # activate to have filtered noise
        if 0:
            
            # filter parameters
            cutoff_freq = 1.0/4.0
            
            ripple_db = 60                      # ripples and transition width of the filter
            width = 1 / 5.0
            
            N_filter, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter

            taps = signal.firwin( N_filter, cutoff=cutoff_freq,  window=('kaiser', beta))
            f_1 = signal.lfilter(taps, 1.0, f_1)
            

        # find periodogram by simple fft and abs()**2
        psd_p[n, :] = own_spec.find_periodogram(f_1, Ome)
        
        psd_bartlett[n, :] = own_spec.find_bartlett_estimate(f_1, M, Ome)
        
        n += 1
    
        # show progress
        done = float(n)/N_trials*100.0
        print 'Done: %3.2f percent' % done
     
    
    # averaging and finding variance
    psd_p_average = psd_p.mean(axis=0)
    #psd_p_average = psd_p_average / np.max(np.abs(psd_p_average))
    
    psd_bartlett_average = psd_bartlett.mean(axis=0)    
    #psd_bartlett_average = psd_bartlett_average / np.max(np.abs(psd_bartlett_average))

    psd_p_deviation = psd_p.var(axis=0)    
    psd_bartlett_deviation = psd_bartlett.var(axis=0)      
     
    #psd_p_deviation = psd_p.std(axis=0) 
    #psd_p_bartlett_deviation = psd_p_bartlett.std(axis=0)     

    

    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    

    # plot psd1
    plt.figure(2)    
    
    plt.subplot(121)    
    plt.plot(Ome, psd_p_average)      
    plt.plot(Ome, psd_p_average - psd_p_deviation)          
    plt.plot(Ome, psd_p_average + psd_p_deviation)          

    plt.title('Periodogram')    
    plt.grid(True); 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_p(\Omega)$')   
    #plt.axis([-2, 2, -1, 7])
    plt.xlim(-1.5, 1.5)    
    plt.ylim( (-6, 20 ) )

    plt.subplot(122)    
    plt.plot(Ome, psd_bartlett_average)  
    plt.plot(Ome, psd_bartlett_average - psd_bartlett_deviation)          
    plt.plot(Ome, psd_bartlett_average + psd_bartlett_deviation)          
    
    plt.title('Bartlett method')
    plt.grid(True); 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_B(\Omega)$')       
    #plt.axis([-2, 2, -1, 7])
    plt.xlim(-1.5, 1.5)    
    plt.ylim( (-6, 20 ) )
    
    plt.show()
        

    


    

    
########################
# make it executable
########################
if __name__ == "__main__":
    main()