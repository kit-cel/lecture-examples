# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 13:25:04 2015

@author: jaekel
"""

########################
######## importing
########################
  
from scipy import signal  
  
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import time


########################
# main function
########################
def main():
    
    # define parameters
    N_trials = 1e2
    N_test = 1e3
    N_step = 5

    performance = np.zeros( (3, N_test) )

    # loop for sequence lengths    
    for _n in np.arange(1,N_test+1, N_step):

        
        #x = np.random.randn(_n)
        #h = np.random.randn(_n)
        
        x = np.ones(_n)
        h = np.ones(_n)
        
        H = np.fft.fft(h)

        # duration of standard convolution        
        tic_1 = time.time()

        for _k in np.arange(N_trials):
            #y = np.convolve(x, h)
            y = signal.lfilter(h, 1, x)

            
        tic_2 = time.time()
        
        for _k in np.arange(N_trials):
            y = np.fft.ifft( np.fft.fft(x) * H)
        
        tic_3 = time.time()
        
        for _k in np.arange(N_trials):
            y = np.fft.ifft( np.fft.fft(x) * np.fft.fft(h) ) 
        
        tic_4 = time.time()
        
        performance[0, _n-1] = (tic_2 - tic_1)/N_trials
        performance[1, _n-1] = (tic_3 - tic_2)/N_trials
        performance[2, _n-1] = (tic_4 - tic_3)/N_trials
        
        print _n
        
    # plotting
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    
   
    
    plt.figure(1)

    plt.plot( np.arange(N_test), performance[0, :], label='conv')
    plt.plot( np.arange(N_test), performance[1, :], label='fast conv., precomp. H')
    plt.plot( np.arange(N_test), performance[2, :], label='fast conv.')    
    
    plt.grid(True)
    plt.legend(loc='upper left')   

    plt.xlabel('Sequence length')        
    plt.ylabel('time $/ \mathrm{s}$')

    plt.show()
    
    # that's it folks
    print('Done!')    
    

   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()        