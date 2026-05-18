# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 16:55:58 2014

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
    N = int( 1e3 )
       
    N_acf = np.arange(-N+1, N, 1)
    
    # number of freq. points and freq. range
    N_freq = 512            
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
    # number of realizations for averaging    
    N_trials = int( 1e2 )
    
    acf_biased = np.empty([N_trials, len(N_acf)])
    acf_unbiased = np.empty([N_trials, len(N_acf)])
    
    psd_p = np.empty([N_trials, N_freq])
    psd_c_biased = np.empty([N_trials, N_freq])
    psd_c_unbiased = np.empty([N_trials, N_freq])
    
    n=0
    while n < N_trials:
    
        # generate noise
        f_1 = np.sqrt(2) * np.random.normal(0.0, 1.0, N)

        # activate to have filtered noise
        if 0:
            
            # filter parameters
            cutoff_freq = 1.0/4.0
            
            ripple_db = 60                      # ripples and transition width of the filter
            width = 1 / 5.0
            
            N_filter, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter

            taps = signal.firwin( N_filter, cutoff=cutoff_freq,  window=('kaiser', beta))
            f_1 = signal.lfilter(taps, 1.0, f_1)
            

        # find acf estimations        
        acf_biased[n, :] = own_spec.est_acf(f_1, 'biased') 
        acf_unbiased[n, :] = own_spec.est_acf(f_1, 'unbiased') 
        
        # find periodogram by simple fft and abs()**2
        psd_p[n, :] = own_spec.find_periodogram(f_1, Ome)
        
        psd_c_biased[n, :] = own_spec.find_correlogram( acf_biased[n,:], Ome)
        psd_c_unbiased[n, :] = own_spec.find_correlogram( acf_unbiased[n,:], Ome)
   
        n += 1
    
        # show progress
        done = float(n)/N_trials*100.0
        print('Done: %3.2f percent' % done)
     
    
    # averaging and finding variance
    acf_biased_average = acf_biased.mean(axis=0)
    acf_unbiased_average = acf_unbiased.mean(axis=0)     
    
    psd_p_average = psd_p.mean(axis=0)
    psd_c_biased_average = psd_c_biased.mean(axis=0)    
    psd_c_unbiased_average = psd_c_unbiased.mean(axis=0)  

    psd_p_deviation = psd_p.var(axis=0)    
    psd_c_biased_deviation = psd_c_biased.var(axis=0)      
    psd_c_unbiased_deviation = psd_c_unbiased.var(axis=0)          
     
    #psd_p_deviation = psd_p.std(axis=0) 
    #psd_c_biased_deviation = psd_c_biased.std(axis=0)     
    #psd_c_unbiased_deviation = psd_c_unbiased.std(axis=0)     

    

    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    

    # plot acfs
    plt.figure(1)    
    plt.subplot(111)
    plt.plot(N_acf, acf_unbiased_average,'b', label='$L=N-k$')      
    plt.plot(N_acf, acf_biased_average,'r', label='$L=N$')
    #plt.plot(N_acf_1, (1-abs(N_acf_1)/N_1)*acf_wgn_unbiased_1,'--')
    plt.xlabel('$k$')
    plt.ylabel('$\hat{r}[k]$')   
    plt.grid(True)    
    plt.title('$N=$'+str(N))
    plt.legend(loc='upper right') 
   
    
    # plot psd1
    plt.figure(2)    
    
    plt.subplot(131)    
    plt.plot(Ome, psd_p_average)      
    plt.plot(Ome, psd_p_average - psd_p_deviation)          
    plt.plot(Ome, psd_p_average + psd_p_deviation)          
    #plt.title('$N=$'+str(N))    
    plt.grid(True); #plt.legend(loc='upper right') 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_p(\Omega)$')   
    #plt.axis([-4, 4, -1, 3])
    
    plt.subplot(132)    
    plt.plot(Ome, psd_c_biased_average, label='$L=N$')  
    plt.plot(Ome, psd_c_biased_average - psd_c_biased_deviation)          
    plt.plot(Ome, psd_c_biased_average + psd_c_biased_deviation)          
    plt.grid(True); plt.legend(loc='upper right') 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')       
    #plt.axis([-4, 4, -1, 3])
    
    plt.subplot(133)        
    plt.plot(Ome, psd_c_unbiased_average, label='$L=N-k$')             
    plt.plot(Ome, psd_c_unbiased_average - psd_c_unbiased_deviation)          
    plt.plot(Ome, psd_c_unbiased_average + psd_c_unbiased_deviation)              
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.grid(True); plt.legend(loc='upper right') 
    #plt.axis([-4, 4, -10, 10])
    
    
    plt.show()
        



    
########################
# make it executable
########################
if __name__ == "__main__":
    main()