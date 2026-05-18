# -*- coding: utf-8 -*-
"""
Created on Tue May 13 09:52:28 2014

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
    
    # sampling time
    t_s = 1.
    
    # designing a filter
    nyq_freq = 1.0 /(t_s*2.0)        # Nyquist and cutoff frequency of the filter

    cutoff_freq = nyq_freq / 2.
    width = nyq_freq / 3.0

    ripple_db = 20                      
    
       
    # find filter and its zeros and poles
    fir = 1
    
    if fir:
        N, beta = signal.kaiserord(ripple_db, width)    # find filter order and beta parameter
        b = signal.firwin( N, cutoff=cutoff_freq,  window=('kaiser', beta), nyq=nyq_freq)
        a = np.append(1, np.zeros(len(b)-1))
    else:
        g_pass = .1
        b, a = signal.iirdesign(cutoff_freq, cutoff_freq+width, g_pass, ripple_db)
    
    
    # get zeros and poles out of the polynomials
    zeros, poles, k = signal.tf2zpk(b, a)  
 
    # decompose zeros and poles
    indices = np.where( abs(zeros) <= 1 )

    zeros_mp_1 = zeros[ indices ]
    zeros_ap = np.delete( zeros, indices ) 

    poles_ap = 1./(np.conjugate(zeros_ap))

    poles_mp = poles
    zeros_mp = np.append( zeros_mp_1, poles_ap)
            
    # frequency responses of the filters
    freqs = np.linspace(-np.pi, np.pi, 512)            

    w, H = my_freq_response(zeros, poles, freqs, H_0 = k)

    w_ap, H_ap = my_freq_response(zeros_ap, poles_ap, freqs)    
    H_0_ap = np.prod( abs(poles_ap) )
    H_ap = H_ap * H_0_ap
    
    w_mp, H_mp = my_freq_response(zeros_mp, poles_mp, freqs, H_0 = k)
    H_mp = H_mp / H_0_ap 
    
    # find numerator and denominator polynomials    
    num, denom = signal.zpk2tf(zeros, poles, k)    
    num_mp, denom_mp = signal.zpk2tf(zeros_mp, poles_mp, k)        

    # find impulse responses    
    n_o, h = signal.dimpulse( [num, denom, k])
    h = np.squeeze(np.array(h))
    h = h / np.linalg.norm(h)
    
    n_o, h_mp = signal.dimpulse( [num_mp, denom_mp, k])
    h_mp = np.squeeze(np.array(h_mp))    
    h_mp = h_mp / np.linalg.norm(h_mp)    
    
    
    # plotting
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    
   
    plt.figure(1)
    
    plt.subplot(331)
    plot_zeros_poles(zeros, poles)  
    plt.title('$H(\cdot)$')          

    plt.subplot(332)
    plot_zeros_poles(zeros_ap, poles_ap)  
    plt.title('$H_{ap}(\cdot)$')              
    
    plt.subplot(333)
    plot_zeros_poles(zeros_mp, poles_mp)          
    plt.title('$H_{mp}(\cdot)$')              
        
    plt.subplot(334)        
    plt.plot(w, abs(H))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H(\Omega)|$')        
    
    plt.subplot(335)        
    plt.plot(w, abs(H_ap))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H_{ap}(\Omega)|$')        
    
    plt.subplot(336)        
    plt.plot(w, abs(H_mp))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H_{mp}(\Omega)|$')        
    
    plt.subplot(337)        
    plt.plot(w, np.angle(H))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H(\Omega)|$')        
    
    plt.subplot(338)        
    plt.plot(w, np.angle(H_ap))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H_{ap}(\Omega)|$')        
    
    plt.subplot(339)        
    plt.plot(w, np.angle(H_mp))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H_{mp}(\Omega)|$')        
    

    
    plt.figure(3)    
    plt.plot(n_o, h, label='non mp')
    plt.plot(n_o, h_mp, label='mp')
    plt.legend()
    plt.grid(True)
    plt.xlabel('$n')
    plt.ylabel('$h[n]$')     
    
    
    plt.show()
    
 


########################
# plotting function for poles and zeros of a discrete-time system
########################
def plot_zeros_poles(zer, pol, axes=[], label_x=[], label_y=[]):
    """
    plots zeros and poles as provided
    
    IN: zeros, poles, axis, label x, label y, grid y/n
    """
   
    z_r=[]
    z_i=[]
    for z in zer:
        z_r.append(z.real)
        z_i.append(z.imag)
        
    plt.plot( z_r, z_i, 'o', markersize=16)
  
    p_r=[]
    p_i=[]         
    for p in pol:
        p_r.append(p.real)
        p_i.append(p.imag)

    plt.plot(p_r, p_i, 'x', markersize=16)
    
    plt.grid(True)
    
    if axes:
        plt.axis(axes)

    if label_x:
        plt.xlabel(label_x)
        
    if label_y:
        plt.ylabel(label_y)        
    
 
    # add unit circle
    u_r=[]
    u_i=[]
    phi = np.linspace(0, 2*np.pi, 2056)
    unit_circle = np.exp(1j*phi)
    for u in unit_circle:
        u_r.append(u.real)
        u_i.append(u.imag)
    plt.plot(u_r,u_i)



########################
# numerically deriving the frequency response out of poles and zeros
# motivated by graphical relation of frequency response and pole-zero diagram
########################    
def my_freq_response(zer, pol, freqs, H_0 = 1):
    """
    gives frequency response of a discrete LTI system when zeros, poles and gain are provided
    
    IN: zeros, poles, H_0, number of freq. points
    
    OUT: H(Omega) at given freqs
    """
    
    if freqs is None:
        num_freq_points = 512
        Omega = np.linspace(-np.pi, np.pi, num_freq_points, endpoint=False)
    elif isinstance(freqs, int):
        num_freq_points = freqs
        Omega = np.linspace(-np.pi, np.pi, num_freq_points, endpoint=False)
    else:
        num_freq_points = len(freqs)
        Omega = freqs
    
    # initialize freq. response
    H = np.zeros(num_freq_points, dtype = complex)
    
    # determine freq. response by distance and angle to zeros and poles
    for k in np.arange(num_freq_points):
        H[k] = np.prod( np.exp(1j*Omega[k])- zer) / np.prod( np.exp(1j*Omega[k])- pol )
    
    return Omega, H_0 * H
   
   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
