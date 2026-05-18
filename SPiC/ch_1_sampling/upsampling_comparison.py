# -*- coding: utf-8 -*-
"""
Created on Mon May  5 09:10:13 2014

@author: jaekel
"""



########################
######## importing
########################
  
from scipy import signal
import time
import numpy as np
import matplotlib.pyplot as plt


########################
# main function
########################
def main():
    
    # define parameters
    t_min=-10.0
    t_max=10.0    
      
    # sampled version 
    t_s = 1.
    
    t_samp = np.arange(t_min, t_max, t_s)

    # determine according Nyquist band and fft
    f_samp = np.arange( -1/(t_s*2.0), 1/(t_s*2.0), 1/(t_max-t_min) )    
    
    # defining the signal
    x_samp = (signal.gausspulse(t_samp, 1, retenv=1))[1] 
    x_samp = x_samp * np.sqrt( 1 / sum(x_samp**2) )
    X_samp = np.fft.fft(x_samp)

    ##########    
    # upsampling parameters
    M = 8
    t_s_up = t_s/(M)
    t_up = np.arange(t_min, t_max, t_s_up)
    f_up = np.arange( -1/(t_s/M*2.0), 1/(t_s/M*2.0), 1/(t_max-t_min) )    


    # number of trials for estimating duration of the algorithms
    numb_trials = 1e3

    # zero padding of fft
    n = 0
    tic = time.time()
    while n < numb_trials:
        X_samp_FD = np.fft.fft(x_samp)
        X_up_FD = np.append( np.append(X_samp_FD[0:len(X_samp_FD)/2], np.zeros((M-1)*len(X_samp_FD))), X_samp_FD[len(X_samp_FD)/2:])
        x_up_FD = np.fft.ifft(X_up_FD)         
        n += 1
    toc = time.time()
    
    print 'Applying zero-padding in the FFT domain: ', toc-tic

    # in the time domain
    # filter design is done before measuring time

    # designing the filter
    nyq_freq = 1.0 /(t_s_up*2.0)        # Nyquist and cutoff frequency of the filter
    cutoff_freq = 1.0/(t_s*2.0)
    ripple_db = 30                      # ripples and transition width of the filter
    width = nyq_freq / 20.0
    N, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter
    taps = signal.firwin( N, cutoff=cutoff_freq,  window=('kaiser', beta), nyq=nyq_freq)

    print '\nNumber of taps in the filtering: ', len(taps)

    n = 0
    tic = time.time()
    while n < numb_trials:
        x_up_TD = np.zeros(( (M)*len(x_samp), ) )
        x_up_TD[::M] = x_samp
        x_up_TD = signal.lfilter(taps, 1.0, x_up_TD)
        n += 1
    toc = time.time()
    
    print 'Applying time domain upsampling and filtering: ', toc-tic

    # fft of TD upsampling for illustration
    X_up_TD = np.fft.fft(x_up_TD)


   # plotting
    if 0:
        plt.close('all')    
        
        plt.figure(1)
        plt.subplot(221)
        plt.stem(t_samp, x_samp)
        plt.stem(t_up, x_up_FD,'r')
        plt.grid(True)
        plt.xlabel('n')
        plt.ylabel('x[n]')
    
        plt.subplot(222)        
        plt.stem(f_samp, np.fft.fftshift(abs(X_samp)**2))
        plt.stem(f_up, np.fft.fftshift(abs(X_up_FD)**2),'r')
        plt.grid(True)
        plt.xlabel('k')
        plt.ylabel('|X[k]|^2')
        
        plt.subplot(223)
        plt.stem(t_samp, x_samp)
        plt.stem(t_up, x_up_TD,'r')
        plt.grid(True)
        plt.xlabel('n')
        plt.ylabel('x[n]')
    
        plt.subplot(224)        
        plt.stem(f_samp, np.fft.fftshift(abs(X_samp)**2))    
        plt.stem(f_up, np.fft.fftshift(abs(X_up_TD)**2),'r')
        plt.grid(True)
        plt.xlabel('k')
        plt.ylabel('|X[k]|^2')
        
        plt.show()
  

########################
# make it executable
########################
if __name__ == "__main__":
    main()
