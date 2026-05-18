# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 11:30:01 2015

@author: jaekel
"""

########################
######## importing
########################
  
from scipy import signal

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    
    # define parameters
    t_min=0.0
    t_max=10.0    
    
    N_fft = 512
      
    # sampled version 
    t_s = .125
    
    t_samp = np.arange(t_min, t_max, t_s)
    x_samp = (signal.gausspulse(t_samp-3, 1, retenv=1))[1] 
    x_samp = x_samp * np.sqrt( 1 / sum(x_samp**2) )

    # determine according Nyquist band and fft
    f_samp = np.arange( -1/(t_s*2.0), 1/(t_s*2.0), 1/(N_fft*t_s) )    
    X_samp = np.fft.fftshift(np.fft.fft(x_samp, N_fft))
    
    # downsampling
    M = 8
    t_s_down = t_s*M
    t_down = np.arange(t_min, t_max, t_s_down)
    
    # upsample signal by factor M using lambda function
    x_down = x_samp[::M]

    # determine according Nyquist band and fft
    f_down = np.arange( -1/(t_s_down*2.0), 1/(t_s_down*2.0), 1/(N_fft*t_s_down) )    
    X_down = np.fft.fftshift(np.fft.fft(x_down, N_fft))
    
    
    # plotting
    font = {'size'   : 36}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    
   
    
    plt.figure(1)
    plt.subplot(221)
    plt.plot(t_samp, x_samp, 'bo', markersize=10, label='$x[n]$')
    plt.grid(True)
    plt.axis(xmin=0, xmax=6, ymin=-.1, ymax=.5)    
    #plt.ylabel('$x[n]$')
    plt.legend(loc='upper right')
            
    plt.subplot(223)
    plt.plot(t_down, x_down, 'bo', markersize=10, label='$x_\\mathrm{down}[n]$')
    plt.grid(True)
    plt.axis(xmin=0, xmax=6, ymin=-.1, ymax=.5)    
    plt.xlabel('$t/\mathrm{s}$')
    #plt.ylabel('$x_\mathrm{down}[n]$')            
    plt.legend(loc='upper right')
          
    plt.subplot(222)        
    plt.plot(f_samp, abs(X_samp)**2, label='$|X[k]|^2$')
    plt.grid(True)
    #plt.ylabel('$|X[k]|^2$')
    #plt.xlim( (-1/(2*t_s), 1/(2*t_s)))
    plt.legend(loc='upper right')
    
    plt.subplot(224)        
    plt.plot(f_down, abs(X_down)**2, label='$|X_\\mathrm{down}[k]|^2$')
    plt.grid(True)
    #plt.xlim( (-1/(2*t_s), 1/(2*t_s)))
    plt.xlabel('$f/\mathrm{Hz}$')
    #plt.ylabel('$|X[k]|^2$')
    plt.legend(loc='upper right')
    
    plt.show()
    
    # that's it folks
    print('Done!')    
    

   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
