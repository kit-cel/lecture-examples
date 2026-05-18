# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 14:13:12 2014

@author: jaekel
"""


########################
# compare different methods for estimation of psd
########################

# importing relevant stuff
import numpy as np
from scipy import signal
import own_utils_est_psd as own_spec

import matplotlib.pyplot as plt
import matplotlib


# choose signal type out of 'noise', 'filtered', 'oscs'
signal_type = 'oscs'


########################
# main function
########################
def main():
    # parameters: number of signal samples
    N = 1e3
    N_vec = np.arange(0, N)
    
    # member per group for Bartlett's method and overlap for Welch
    M = N/10.0
    O = M/2.0
       
    # number of freq. points and freq. range
    N_freq = 512            
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
    # number of realizations for averaging    
    N_trials = 1e2
    
    # initializa arrays    
    psd_p = np.empty([N_trials, N_freq])
    psd_blackman_tukey = np.empty([N_trials, N_freq])        
    psd_bartlett = np.empty([N_trials, N_freq])            
    psd_welch = np.empty([N_trials, N_freq])    
    
    hamming = np.append( np.zeros(N-M), np.append(signal.blackman(2*M), np.zeros( N-1-M)))    

    # loop for realizations, averaging and variance
    n=0
    while n < N_trials:

        # generate signal    
        if signal_type == 'noise':        

            f_1 = np.sqrt(2) * np.random.normal(0.0, 1.0, N)
        
        elif signal_type == 'oscs':
                Omega_0 = 1.0; Omega_1 = 1.2
                f_1 = np.sin(Omega_0 * N_vec) + np.sin(Omega_1 * N_vec)  + np.random.normal(0.0, 1.0, N)
    
        elif signal_type == 'filtered':
                f_1 = np.sqrt(2) * np.random.normal(0.0, 1.0, N)
                
                # filter parameters
                cutoff_freq = 1.05/np.pi
                
                ripple_db = 60                      # ripples and transition width of the filter
                width = 1 / 10.0
                
                N_filter, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter
    
                taps = signal.firwin( N_filter, cutoff=cutoff_freq,  window=('kaiser', beta))
                f_1 = signal.lfilter(taps, 1.0, f_1)
            

        # find different estimators
        psd_p[n, :] = own_spec.find_periodogram(f_1, Ome)
        
        acf = own_spec.est_acf(f_1, 'biased') 
        psd_blackman_tukey[n, :] = own_spec.find_correlogram(acf * hamming, Ome)
        
        psd_bartlett[n, :] = own_spec.find_bartlett_estimate(f_1, M, Ome)        
        
        psd_welch[n, :] = own_spec.find_welch_estimate(f_1, M, O, Ome)
        
        #psd_welch_onboard[n, :] = np.fft.fftshift((signal.welch(f_1, window='blackman', nfft= N_freq, nperseg=M, noverlap=O/2, return_onesided=False))[1])
        
        n += 1
    
        # show progress
        done = float(n)/N_trials*100.0
        print 'Done: %3.2f percent' % done
     
    
    # averaging and finding variance
    psd_p_average = psd_p.mean(axis=0)    
    psd_blackman_tukey_average = psd_blackman_tukey.mean(axis=0)        
    psd_bartlett_average = psd_bartlett.mean(axis=0)        
    psd_welch_average = psd_welch.mean(axis=0)    

    if signal_type == 'filtered' or signal_type == 'oscs':
        psd_p_average = psd_p_average / np.max(psd_p_average)
        psd_blackman_tukey_average = psd_blackman_tukey_average / np.max(psd_blackman_tukey_average)
        psd_bartlett_average = psd_bartlett_average / np.max(psd_bartlett_average )
        psd_welch_average = psd_welch_average / np.max(psd_welch_average) 
    

    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 30}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    

    # plot psd1
    plt.figure(13)    
    
    plt.subplot(221)
    
    if signal_type == 'oscs':# or signal_type == 'filtered':
        plt.plot(Ome, 10*np.log10(psd_p_average), label='Period.')    
        plt.axis([-4, 4, -30,0])
        plt.ylabel('$\Phi(\Omega)$ (dB)')
        
    else:
        plt.plot(Ome, psd_p_average, label='Period.')      
        plt.axis([-4, 4, 0, 10])
   
        plt.ylabel('$\Phi(\Omega)$')
    plt.grid(True)
    plt.title('Period.')
    
    
    plt.subplot(222)    
    if signal_type == 'oscs':# or signal_type == 'filtered':
        plt.plot(Ome, 10*np.log10(psd_blackman_tukey_average,), label='Blackman-T.')    
        plt.axis([-4, 4, -30,0])
        plt.ylabel('$\Phi(\Omega)$ (dB)')
        
    else:
        plt.plot(Ome, psd_blackman_tukey_average, label='Blackman-T..')      
        plt.axis([-4, 4, 0, 10])
         
    plt.grid(True); 
    plt.title('Blackman-T.')

    

    plt.subplot(223)
    if signal_type == 'oscs':# or signal_type == 'filtered':
        plt.plot(Ome, 10*np.log10(psd_bartlett_average,), label='Bartlett')    
        plt.axis([-4, 4, -30,0])
        plt.ylabel('$\Phi(\Omega)$ (dB)')
        
    else:
        plt.plot(Ome, psd_bartlett_average, label='Bartlett.')      
        plt.axis([-4, 4, 0, 10])
        
    plt.grid(True); 
    plt.title('Bartlett')    
    plt.xlabel('$\Omega$')
    plt.ylabel('$\Phi(\Omega)$')
    
    plt.subplot(224)
    if signal_type == 'oscs':# or signal_type == 'filtered':
        plt.plot(Ome, 10*np.log10(psd_welch_average), label='Welch')    
        plt.axis([-4, 4, -30,0])
        plt.ylabel('$\Phi(\Omega)$ (dB)')
        
    else:
        plt.plot(Ome, psd_welch_average, label='Welch')      
        plt.axis([-4, 4, 0, 10])
        
        
    plt.grid(True); 
    plt.title('Welch')
    plt.xlabel('$\Omega$')
    #plt.axis([-4, 4, 0, 20])
    #plt.legend(loc=1)

    plt.show()
        

    
########################
# make it executable
########################
if __name__ == "__main__":
    main()