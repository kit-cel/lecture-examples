# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 15:46:00 2015

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
    K_1 = 21
    K_2 = 105    
    
    # define ideal lowpass in the frequency regime
    H_w = np.zeros(N_fft)
    H_w[ np.where( np.abs(f)<f_g) ] = 1
    
    # find impulse response by IFFT and restricting to K values
    h_1_part = np.fft.ifft( H_w*np.exp(-1j*2*np.pi*f*(K_1-1)/2), N_fft)[:(K_1+1)//2]
    h_1 = np.append( h_1_part, (h_1_part[::-1])[1:])
    #h /= np.linalg.norm(h)
    
    h_2_part = np.fft.ifft( H_w*np.exp(-1j*2*np.pi*f*(K_2-1)/2), N_fft)[:(K_2+1)//2]
    h_2 = np.append( h_2_part, (h_2_part[::-1])[1:])#, h_2_part)
    #h_2 /= np.linalg.norm(h_2)        
    
    
    # find frequency responses
    #H = np.fft.fftshift( find_magnitude(h, 2*np.pi*f*t_s) ) 
    #H /= H[N_fft/2]
    
    #H_2 = np.fft.fftshift( find_magnitude(h_2, 2*np.pi*f*t_s) ) 
    #H_2 /= H_2[N_fft/2]

    freq, H_1 = signal.freqz(h_1, worN=f*2*np.pi, whole=True) 
    H_1 = np.fft.fftshift(H_1)
    
    freq, H_2 = signal.freqz(h_2, worN=f*2*np.pi, whole=True)
    H_2 = np.fft.fftshift(H_2)
    
    # plotting
    plt.figure(1)
    plt.clf()
    
    plt.subplot(131)
    plt.plot( np.arange(len(h_1)), h_1, label='$K=$'+str(K_1))
    plt.plot( np.arange(len(h_2)), h_2, label='$K=$'+str(K_2))    
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.xlabel('$n$')
    plt.title('$h[n]$')
  
    plt.subplot(132)
    plt.plot( f, np.abs(H_1), label='$K=$'+str(K_1))
    plt.plot( f, np.abs(H_2), label='$K=$'+str(K_2))    
    plt.grid(True)   
    plt.legend(loc='center right')    
    plt.xlabel('$f/\mathrm{Hz}$')
    plt.title('$|H(f)|$')    

    plt.subplot(133)
    plt.plot( f, 10*np.log10(np.abs(H_1)), label='$K=$'+str(K_1))
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
