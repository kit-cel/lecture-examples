# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:38:26 2016

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
    
    N_fft = 5096
    f = np.arange( -f_s/2, f_s/2, f_s/N_fft)    

    # end of passband and filter length
    f_g = f_s/3
    K = 105


    #
    # design Fourier approximation
    #
    
    # define ideal lowpass in the frequency regime
    H_w = np.zeros(N_fft)
    H_w[ np.where( np.abs(f)<f_g) ] = 1
    
    # find impulse response by IFFT and restricting to K values
    h_F_part = np.fft.ifft( H_w * np.exp( -1j * 2 * np.pi * f * (K-1)/2), N_fft)[:(K+1)/2]
    h_F = np.append( h_F_part, ( h_F_part[::-1] )[1:] )
    #h /= np.linalg.norm(h)  
    
    # find frequency responses
    freq, H_F = signal.freqz( h_F, worN=f*2*np.pi, whole=True ) 
    H_F = np.fft.fftshift( H_F )
    
    #
    # design Remez filter
    #    
    
    # find impulse response by IFFT and restricting to K values
    h_R = signal.remez( K, [0, .9*f_g, 1.1*f_g, .5], [1,0], weight=[1,1] )
    
    # find frequency responses
    freq, H_R = signal.freqz( h_R, worN=f*2*np.pi, whole=True )  
    #H = np.fft.fftshift(H)
    
    
    #
    # determine errors
    #

    print('Fourier approx. with K={} gives quad. error: {:2.2f}\n'.format( K, np.linalg.norm( H_w - H_F )**2 ) )
    print('Remez approx. with K={} gives quad. error: {:2.2f}\n'.format( K, np.linalg.norm( H_w - H_R )**2 ) )
    
    #print(h_F, h_R)
#    plt.figure(3)
#    plt.plot( f, np.abs( np.real( H_w-H_F )), label='F.' )
#    plt.plot( f, np.abs( np.real( H_w-H_R )), label='R.' )
#    plt.grid( True )
#    plt.legend( loc='upper right')
    
    
    # plotting
    plt.figure(1)
    plt.clf()
    
    plt.subplot(131)
    plt.plot( np.arange( len(h_F) ), np.real( h_F ), label='F., $K=$'+str(K))
    plt.plot( np.arange( len(h_R) ), np.real( h_R ), label='R., $K=$'+str(K))
      
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.xlabel('$n$')
    plt.title('$h[n]$')
  
  
    plt.subplot(132)
    plt.plot( f, np.abs( H_F), label='F., $K=$'+str(K) )
    plt.plot( f, np.abs( H_R), label='R., $K=$'+str(K) )
    plt.plot( f, np.abs( H_w), label='W., $K=$'+str(K) )    
    
    plt.grid(True)   
    plt.legend(loc='lower right')
    plt.xlabel('$f/\mathrm{Hz}$')
    plt.title('$|H(f)|$')    
   
    plt.subplot(133)
    plt.plot( f, 10*np.log10( np.abs( H_F) ), label='F., $K=$'+str(K) )
    plt.plot( f, 10*np.log10( np.abs( H_R) ), label='R., $K=$'+str(K) )
    plt.plot( f, 10*np.log10( np.abs( H_w) ), label='W., $K=$'+str(K) )    
    
    plt.grid(True)   
    plt.legend(loc='lower right')
    plt.xlabel('$f/\mathrm{Hz}$')
    plt.title('$|H(f)|$ (dB)') 
    
    
#    plt.subplot(133)
#    plt.plot( f, 10*np.log10( np.abs( H_F )**2 ), label='F., $K=$'+str(K) )   
#    plt.plot( f, 10*np.log10( np.abs( H_R )**2 ), label='R., $K=$'+str(K) )
#    plt.plot( f, 10*np.log10( np.abs( H_w )**2 ), label='W., $K=$'+str(K) )
#
#    plt.grid(True)   
#    plt.legend(loc='lower right')   
#    plt.xlabel('$f/\mathrm{Hz}$')
#    plt.title('$|H(f)|^2 \\; (dB)$') 
   

   
    plt.show()
    
    # that's it folks
    print('Done!')

    
    
########################
# make it executable
########################
if __name__ == "__main__":
    main()
