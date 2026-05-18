# -*- coding: utf-8 -*-
"""
Created on Thu Apr 09 13:53:43 2015

@author: jaekel
"""


from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from scipy import signal

# nice tex plots
font = {'size'   : 26}
plt.rc('font', **font)
plt.rc('text', usetex=True)


 
########################
# main function
########################
def main():
    
    # sampling time
    t_s = 1.
    f_s = 1. / t_s
    
    N_fft = 2048
    f = np.arange( -f_s/2, f_s/2, f_s/N_fft)    

    # end of passband and filter length
    f_g = f_s/3
    K_1 = 11
    K_2 = 51    
    
    # define ideal lowpass in the frequency regime
    H_w = np.zeros(N_fft)
    H_w[ np.where( np.abs(f)<f_g) ] = 1
    
    # find impulse response 
    h = signal.remez( K_1, [0, f_g, 1.2*f_g, .5], [1,0])
    h_2 = signal.remez(K_2, [0, f_g, 1.2*f_g, .5], [1,0])    
    
    
    # find frequency responses
    freq, H = signal.freqz(h, worN=f*2*np.pi, whole=True) 
    #H = np.fft.fftshift(H)
    
    freq, H_2 = signal.freqz(h_2, worN=f*2*np.pi, whole=True)
    #H_2 = np.fft.fftshift(H_2)
    
    # plotting
    plt.figure(1)
    plt.clf()
    
    plt.subplot(131)
    plt.plot( np.arange(len(h)), h, label='$K=$'+str(K_1))
    plt.plot( np.arange(len(h_2)), h_2, label='$K=$'+str(K_2))    
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.xlabel('$n$')
    plt.title('$h[n]$')
  
    plt.subplot(132)
    plt.plot( f, np.abs(H), label='$K=$'+str(K_1))
    plt.plot( f, np.abs(H_2), label='$K=$'+str(K_2))    
    plt.grid(True)   
    plt.legend(loc='center right')    
    plt.xlabel('$f/\mathrm{Hz}$')
    plt.title('$|H(f)|$')    
   
    plt.subplot(133)
    plt.plot( f, 10*np.log10(np.abs(H)), label='$K=$'+str(K_1))
    plt.plot( f, 10*np.log10(np.abs(H_2)), label='$K=$'+str(K_2))    
    plt.grid(True)   
    plt.legend(loc='center right')    
    plt.xlabel('$f/\mathrm{Hz}$')
    plt.title('$|H(f)| \\; (dB)$') 
   
    plt.show()
    
    # that's it folks
    print('Done!')

    
    
    
    
########################
# make it executable
########################
if __name__ == "__main__":
    main()
