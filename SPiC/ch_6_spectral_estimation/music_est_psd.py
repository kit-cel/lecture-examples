# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 12:57:11 2014

@author: jaekel
"""


########################
# find MUSIC estimate
########################

# importing relevant stuff
import numpy as np
import own_utils_est_psd as own_spec 
from scipy import signal

import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    
    # parameters: number of signal samples
    N = int( 1e3 )
    N_vec = np.arange(0, N)
    
    # parameter (system order) used in the MUSIC algorithm
    # number of observations used for the algorithm -->
    # large m dramatically reduce performance
    n = 50
    m = N   
    
    # number of freq. points and freq. range
    N_freq = 512
    Ome = np.linspace(-np.pi, np.pi, N_freq)

    # number of realizations for averaging    
    N_trials = int( 1e2 )
    
    # initializa arrays    
    psd_p = np.empty([N_trials, N_freq])
    psd_music = np.empty([N_trials, N_freq])        

    # frequencies of sinusoids with or w/o leakage
    Omega_0 = 1.0 #Ome[np.random.randint(0, N_freq)]
    Omega_1 = 1.2 #Ome[np.random.randint(0, N_freq)]
    sigma2 = 2

    # loop for realizations, averaging and variance
    k=0
    while k < N_trials:
  
        # different signals to be analyzed      
        choice = 6
        #     
        if choice==1: # AR spectrum out of the "spectrum" homepage
            a = [1, -2.2137, 2.9403, -2.1697, 0.9606]
            f_1 = signal.lfilter([1], a, np.random.normal(0.0, 1.0, N))
            
        elif choice==2: # just noise
            f_1 = np.sqrt(sigma2) * np.random.normal(0.0, 1.0, N)
            
        elif choice==3: # AR spectrum out of Kroschel
            b = [1]
            a = [1 -1.372, 1.843, -1.238, .849]
            f_1 = signal.lfilter(b, a, np.sqrt(sigma2)*N_vec)     
            
        elif choice==4: # broad spectrum out of Kroschel
            b = [1, 0, 0, 0, -.5]
            a = [1]
            f_1 = signal.lfilter(b, a, np.sqrt(sigma2)*N_vec)
            
        elif choice==5: # two sinusoids with noise
            f_1 = np.sin(Omega_0 * N_vec) + np.sin(Omega_1 * N_vec)  + np.sqrt(sigma2)*np.random.normal(0.0, 1.0, N)
      
        elif choice==6: # two complex sinusoids with noise
            f_1 = np.exp(1j*Omega_0 * N_vec) + np.exp(1j*Omega_1 * N_vec)  + np.sqrt(sigma2)*np.random.normal(0.0, 1.0, N)
        
        # find estimators
        psd_p[k, :] = own_spec.find_periodogram(f_1, Ome) #
        psd_music[k, :] = own_spec.find_music_estimate(f_1, n, N, Ome)
        
        k += 1
    
        # show progress
        done = float(k)/N_trials*100.0
        print('Done: %3.2f percent' % done)
     
    
    # averaging and finding variance
    psd_p_average = psd_p.mean(axis=0)
    psd_p_average = psd_p_average/np.max(psd_p_average)
    
    psd_music_average = psd_music.mean(axis=0)
    psd_music_average = psd_music_average/np.max(psd_music_average)        
    
    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    

    # plot psd1
    plt.figure(13)    
    
    plt.subplot(121)
    #plt.plot(Ome, psd_p_average)              
    #plt.plot(Ome, psd_p_average, 'x')          
    plt.plot(Ome, 10*np.log10(psd_p_average))      
    plt.ylabel('$\Phi(\Omega)$ (dB)')
    plt.grid(True)
    plt.title('Period.')
    plt.axis([-4, 4, -30, 0])
    plt.xlabel('$\Omega$')    
    
    plt.subplot(122)    
    #plt.plot(Ome, psd_music_average)            
    #plt.plot(Ome, psd_music_average, 'x')        
    plt.plot(Ome,  10*np.log10(psd_music_average))        
    plt.grid(True) 
    plt.title('MUSIC; order = '+str(n))
    plt.axis([-4, 4, -30, 0])
    #plt.xlim(0, np.pi)
    plt.xlabel('$\Omega$')

    plt.show()
        
    
    
########################
# make it executable
########################
if __name__ == "__main__":
    main()