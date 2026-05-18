# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 08:49:43 2015

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
    t_min=-10.0
    t_max=10.0    
      
    # sampled version 
    t_s = 1.0
    
    t_samp = np.arange(t_min, t_max, t_s)
    
    x_samp = (signal.gausspulse(t_samp/3, 1, retenv=1))[1] 
    x_samp = x_samp * np.sqrt( 1 / sum(x_samp**2) )

    # define several impulse responses
    K = 10
    length = 20

    h_1 = np.zeros(length)
    h_1[K] = 1
    
    h_2 = np.zeros(length)
    h_2[:K] = 1
    h_2 /= np.linalg.norm(h_2)
   
    # get outputs
    y_1 = np.convolve(x_samp, h_1)
    y_2 = np.convolve(x_samp, h_2)    
    
       
    # spectra
    X = np.fft.fftshift( np.fft.fft( x_samp, 512))
    Y_1 = np.fft.fftshift( np.fft.fft( y_1, 512))
    Y_2 = np.fft.fftshift( np.fft.fft( y_2, 512))
    
    f = np.arange(-1/(2*t_s), 1/(2*t_s), 1/(t_s*512))
    
    # plotting
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    

    
    plt.figure(1)
    plt.clf()

    plt.subplot(231)
    plt.plot( np.arange(len(h_1)), h_1, 'o', label='$h[n]$')
    plt.grid(True)
    plt.axis([0, length, -.1, 1.1])
    plt.legend( loc = 'upper right')

    plt.subplot(232)
    plt.plot( np.arange(len(x_samp)), x_samp, 'o', label='$x[n]$')    
    plt.plot( np.arange(len(y_1)), y_1, 'o', label='$y[n ]$')
    plt.grid(True)
    plt.legend( loc = 'upper right')

    plt.subplot(233)
    plt.plot(f, np.abs(X), label='$|X(f)|$')     
    plt.plot(f, np.abs(Y_1), label='$|Y(f)|$') 
    plt.grid(True)
    plt.legend(loc='upper right')        
    #plt.title('$|X(f)|, |Y(f)|$')

    plt.subplot(234)
    plt.plot( np.arange(len(h_2)), h_2, 'o', label='$h[n]$')
    plt.grid(True)
    plt.axis([0, length, -.1, 1.1])
    plt.xlabel('$n$')     
    plt.legend( loc = 'upper right')
      
    plt.subplot(235)
    plt.plot( np.arange(len(x_samp)), x_samp, 'o', label='$x[n]$')    
    plt.plot( np.arange(len(y_2)), y_2, 'o', label='$y[n ]$')
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.xlabel('$n$') 
    
    plt.subplot(236)
    plt.plot(f, np.abs(X), label='$|X(f)|$')     
    plt.plot(f, np.abs(Y_2), label='$|Y(f)|$') 
    plt.grid(True)
    plt.legend(loc='upper right')    
    plt.xlabel('$f/\mathrm{Hz}$')    

  
    
     
   
    plt.show()
    
    # that's it folks
    print('Done!')
   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
