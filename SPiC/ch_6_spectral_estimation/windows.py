# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 14:42:29 2014

@author: jaekel
"""

########################
# illustrating asymptotically unbiased estimation of psd
########################

# importing relevant stuff
import numpy as np
from scipy import signal

import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    # parameters: number of points in first resp. second function and in frequency domain
    N = 120
    M = 51
        
    n = np.arange(-N/2, N/2)
    
    
    start = int( N/2-(M-1)/2 )
    
    # define windows
    rect = np.zeros(N)
    rect[ start : start+M ] = 1
    
    win_tria = signal.triang(M)
    #tria = np.concatenate( (np.zeros( N/2-(M-1)/2), win_tria, np.zeros( N/2-(M-1)/2-1) ) )
    tria = np.zeros(N)
    tria[ start:start+M] = win_tria
    
    win_hann = signal.hann(M)
    hann = np.zeros(N)
    hann[ start: start+M] = win_hann

    win_hamming = signal.hamming(M)
    hamming = np.zeros(N)
    hamming[ start:start+M] = win_hamming
    
    win_blackman = signal.blackman(M)
    blackman = np.zeros(N)
    blackman[ start: start+M] = win_blackman


    # frequency range, applying zero-padding
    zp = 8
    
    Ome = np.linspace(-np.pi, np.pi, N*zp)
    
    RECT = find_periodogram( np.append(rect, np.zeros((zp-1)*len(n))), Ome)
    RECT = RECT / np.max(RECT) 
    
    TRIA = find_periodogram( np.append(tria, np.zeros((zp-1)*len(n))), Ome)
    TRIA = TRIA / np.max(TRIA)
    
    HANN = find_periodogram( np.append(hann, np.zeros((zp-1)*len(n))), Ome)
    HANN = HANN / np.max(HANN)
    
    HAMMING = find_periodogram( np.append(hamming, np.zeros((zp-1)*len(n))), Ome)
    HAMMING = HAMMING / np.max(HAMMING)
    
    BLACKMAN = find_periodogram( np.append(blackman, np.zeros((zp-1)*len(n))), Ome)
    BLACKMAN = BLACKMAN / np.max(BLACKMAN)    

    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    

    # plot windows
    plt.figure(2)    
    
    plt.subplot(221)   

    plt.plot(n, rect, linewidth=2.0, label='rect.')
    plt.plot(n, tria, linewidth=2.0, label='triang.')
    plt.plot(n, hann, linewidth=2.0, label='Hann')

    #plt.xlabel('$k$')
    plt.ylabel('$w[k]$')   
    plt.grid(True)    
    plt.legend(loc='upper right')      
    plt.axis([-N/2+1, N/2+1, 0, 1.2]) 
      
    plt.subplot(223)   
    plt.plot(n, rect, linewidth=2.0, label='rect.')
    plt.plot(n, hamming, linewidth=2.0, label='Hamming')    
    plt.plot(n, blackman, linewidth=2.0, label='Blackman')
    
    plt.xlabel('$k$')
    plt.ylabel('$w[k]$')   
    plt.grid(True)    
    plt.legend(loc='upper right')  
    plt.axis([-N/2+1, N/2+1, 0, 1.2]) 
 
    plt.subplot(222)   

    plt.plot(Ome, 10*np.log10(RECT), linewidth=2.0, label='rect.')
    plt.plot(Ome, 10*np.log10(TRIA), linewidth=2.0, label='triang.')
    plt.plot(Ome, 10*np.log10(HANN), linewidth=2.0, label='Hann')    

    plt.ylabel('$|W(\Omega)|^2$ (dB)')   
    plt.grid(True)    
    plt.legend(loc='upper right')    
    plt.axis([-1, 1, -80, 10])     
    
    plt.subplot(224)
    plt.plot(Ome, 10*np.log10(RECT), linewidth=2.0, label='rect.')
    plt.plot(Ome, 10*np.log10(HAMMING), linewidth=2.0, label='Hamming')     
    plt.plot(Ome, 10*np.log10(BLACKMAN), linewidth=2.0, label='Blackman')     
    
    plt.xlabel('$\Omega$')
    plt.ylabel('$|W(\Omega)|^2$ (dB)')   
    plt.grid(True)    
    plt.legend(loc='upper right')    
    plt.axis([-1, 1, -80, 10]) 


    plt.show()
        


########################
# own periodogram estimator
########################
def find_periodogram(y, omega):
    """
    estimates periodogram out of the given observation at the frequencies specified in omega
    
    IN: observation y, frequencies
    OUT: psd estimator
    """

    
    #periodogram = 1./len(y) * abs(  np.fft.fft( y) )**2
    #periodogram = np.fft.fftshift( periodogram )

    N = len(y)

    per = (0+0*1j)*omega  #np.array( omega, dtype=complex)
        
    for p in np.arange(0, N):
        per += y[p]*np.exp(-1j*omega*(p+1))
        
    per = (abs(per)**2)/N
        
    return per


    
########################
# make it executable
########################
if __name__ == "__main__":
    main()