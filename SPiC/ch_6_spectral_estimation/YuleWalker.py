# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 13:37:51 2014

@author: jaekel
"""


########################
# illustrating Yule-Walker method
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
    # parameters: number of signal samples
    N = int( 1e3 )
    N_vec = np.arange(0, N)
    
    # order of Yule-Walker approach; q<N required
    q = 3
    q_2 = 50
       
    # variance of noise used in the illustration
    sigma2 = 1
       
    # number of freq. points and freq. range
    N_freq = 512            
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
    # number of realizations for averaging    
    N_trials = 1
    
    # initializa arrays    
    #psd_p = np.empty([N_trials, N_freq])
    psd_p = np.empty([N_trials, N_freq])    
    psd_yulewalker = np.empty([N_trials, N_freq])        
    psd_yulewalker_2 = np.empty([N_trials, N_freq])   

    Omega_0 = 1.0
    Omega_1 = 1.2
    sigma2 = 2

    # loop for realizations, averaging and variance
    n=0
    while n < N_trials:
  
        # different signals to be analyzed      
        choice = 6
        #     
        if choice==1: # AR spectrum out of the "spectrum" homepage
            a = [1, -2.2137, 2.9403, -2.1697, 0.9606]
            f_1 = signal.lfilter([1], a, np.random.normal(0.0, 1.0, N))
        elif choice==2: # just noise
            f_1 = np.sqrt(sigma2) * np.random.normal(0.0, 1.0, N)
        elif choice==3: # AR spectrum out of Kroschel
            b = [1]; a = [1 -1.372, 1.843, -1.238, .849]; sigma2 = .0032
            f_1 = signal.lfilter(b, a, np.sqrt(sigma2)*Omega_1*N_vec)            
        elif choice==4: # broad spectrum out of Kroschel
            b = [1, 0, 0, 0, -.5]; a = [1]; sigma2 = .44
            f_1 = signal.lfilter(b, a, np.sqrt(sigma2)*Omega_1*N_vec)
        elif choice==5: # two sinusoids with noise
            Omega_0 = 1.0; Omega_1 = 1.2; sigma2 = 2
            f_1 = np.sin(Omega_0 * N_vec) + np.sin(Omega_1 * N_vec)  + np.sqrt(sigma2)*np.random.normal(0.0, 1.0, N)
        elif choice==6: # two complex sinusoids with noise
            f_1 = np.exp(1j*Omega_0 * N_vec) + np.exp(1j*Omega_1 * N_vec)  + np.sqrt(sigma2)*np.random.normal(0.0, 1.0, N)
            
        # activate to filter the signal
        if 0:
            # filter parameters
            cutoff_freq = 1.0/np.pi
            ripple_db = 30                      # ripples and transition width of the filter
            width = 1 / (5.0 * np.pi)
            
            N_filter, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter
            taps = signal.firwin( N_filter, cutoff=cutoff_freq,  window=('kaiser', beta))
            f_1 = signal.lfilter(taps, 1.0, f_1)
            

        # find estimators
        psd_p[n, :] = own_spec.find_periodogram(f_1, Ome) #

        # YW of 2 different orders
        psd_yulewalker[n, :] = own_spec.find_yulewalker(f_1, q, sigma2, Ome)    
        psd_yulewalker_2[n, :] = own_spec.find_yulewalker(f_1, q_2, sigma2, Ome)    
        
        # compare to results of lib "spectrum" in order to verify own results
        #AR, P, k = spec.aryule(f_1, 4)
        #psd_yulewalker_2[n, :] = spec.arma2psd(AR, NFFT=N_freq)          

        n += 1
    
        # show progress
        done = float(n)/N_trials*100.0
        print('Done: %3.2f percent' % done)
     
    
    # averaging and finding variance
    psd_p_average = psd_p.mean(axis=0)
    psd_p_average = psd_p_average/np.max(psd_p_average)
    
    psd_yulewalker_average = psd_yulewalker.mean(axis=0)    
    psd_yulewalker_average = psd_yulewalker_average/np.max(psd_yulewalker_average)    
    
    psd_yulewalker_2_average = psd_yulewalker_2.mean(axis=0)        
    psd_yulewalker_2_average = psd_yulewalker_2_average/np.max(psd_yulewalker_2_average)    

    

    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    

    # plot psd1
    plt.figure(13)    
    
    plt.subplot(131)
    #plt.plot(Ome, psd_p_average, label='Period.')          
    plt.plot(Ome, 10*np.log10(psd_p_average))      
    plt.ylabel('$\Phi(\Omega)$ (dB)')
    plt.grid(True)
    plt.title('Periodogram')
    #plt.axis([-4, 4, 0, 20])
    plt.xlabel('$\Omega$')    
    
    plt.subplot(132)    
    #plt.plot(Ome, psd_yulewalker_average, label='Yule-Walker')        
    plt.plot(Ome,  10*np.log10(psd_yulewalker_average))        
    plt.grid(True) 
    plt.title('Yule-Walker; order = '+str(q))
    #plt.axis([-4, 4, 0, 20])
    #plt.xlim(0, np.pi)
    plt.xlabel('$\Omega$')
    
    plt.subplot(133)    
    #plt.plot(Ome, psd_yulewalker_2_average, label='Yule-Walker')        
    plt.plot(Ome,  10*np.log10( psd_yulewalker_2_average))
    plt.grid(True) 
    plt.xlabel('$\Omega$')    
    plt.title('Yule-Walker; order = '+str(q_2))

    plt.show()
        


 
    
    
 



   


    

    
########################
# make it executable
########################
if __name__ == "__main__":
    main()