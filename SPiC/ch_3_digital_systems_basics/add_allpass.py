# -*- coding: utf-8 -*-
"""
Created on Mon May 12 16:49:09 2014

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
    
    # first system being minimum phase
    zeros = np.array([0, 0, 0])
    poles = np.array([.5, .25+.8j, .25-.8j])
    num, denom = signal.zpk2tf(zeros, poles, 1)
    
    n_o, h = signal.dimpulse( [num, denom, 1])
    h = np.squeeze(np.array(h))
    
    # freq. range and freq. response
    freqs = np.linspace(-np.pi, np.pi, 512)    
    w, H = signal.freqz(num, denom, worN=freqs, whole=True)    
    H = np.squeeze(np.array(H))    
    
    
    # add all-pass part and find all parameters of this second system
    poles_add = np.array([-.3+.4j, -.3-.4j, -.8, -.1+.7j, -.1-.7j], dtype=complex)
    zeros_add = 1./poles_add.conj()
    
    zeros_2 = np.append(zeros, zeros_add)    
    poles_2 = np.append(poles, poles_add)
    num_2, denom_2 = signal.zpk2tf(zeros_2, poles_2, 1)    
    
    n_o_2, h_2 = signal.dimpulse( [num_2, denom_2, 1])
    h_2 = np.squeeze(np.array(h_2))

    # freq. response
    w_2, H_2 = signal.freqz(num_2, denom_2, worN=freqs, whole=True)    
    H_2 = np.squeeze(np.array(H_2)) 
    
    # partial energy of the 2 systems
    E_h = np.cumsum( abs(h)**2/np.sum(abs(h)**2))
    E_h_2 = np.cumsum( abs(h_2)**2/np.sum(abs(h_2)**2))    
    
    
    # my own freq. reponse
    w_own, H_own = my_freq_response(zeros_2, poles_2, freqs)
        
    
    # plotting
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)  

    plt.figure(1)
    
    plt.subplot(221)
    plot_zeros_poles(zeros, poles, [-1.1, 1.1, -1.1, 1.1], 'Re$\{z\}$', 'Im$\{z\}$')            
            
    plt.subplot(222)
    plt.stem(n_o, h.real)#, n_o, h.imag, 'o')
    #plt.axis([-.1, 30.1, -.1, 1.1])
    plt.grid(True)        
    plt.xlabel('$n$')
    plt.ylabel('$h_{mp}[n]$')
        
    plt.subplot(223)
    plot_zeros_poles(zeros_2, poles_2, [-1.1, 1.1, -1.1, 1.1], 'Re$\{z\}$', 'Im$\{z\}$')             
            
    plt.subplot(224)
    plt.stem(n_o_2, h_2.real)#, n_o, h.imag, 'o')
    #plt.axis([-.1, 30.1, -.1, 1.1])
    plt.grid(True)        
    plt.xlabel('$n$')
    plt.ylabel('$h_2[n]$')    
    
    
    
    plt.figure(2)    
    plt.subplot(321)
    plt.plot(freqs, abs(H))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    #plt.xlabel('$\Omega$')
    plt.ylabel('$|H_{mp}(\Omega)|$')

    plt.subplot(322)    
    plt.plot(freqs, np.angle(H))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    #plt.xlabel('$\Omega$')
    plt.ylabel('$\phi_{mp}(\Omega)$')    
    

    plt.subplot(323)
    plt.plot(freqs, abs(H_2))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    #plt.xlabel('$\Omega$')
    plt.ylabel('$|H_2(\Omega)|$')

    plt.subplot(324)    
    plt.plot(freqs, np.angle(H_2))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    #plt.xlabel('$\Omega$')
    plt.ylabel('$\phi_2(\Omega)$') 

    plt.subplot(325)
    plt.plot(freqs, abs(H_own))
    plt.grid(True)
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H_{2, own}(\Omega)|$') 
    
    plt.subplot(326)    
    plt.plot(freqs, np.angle(H_own))
    plt.grid(True)
    plt.xlabel('$\Omega$')
    plt.ylabel('$\phi_{2, own}(\Omega)$')    
                
    plt.figure(3)    
    plt.plot(n_o, E_h, label='mp')
    plt.plot(n_o, E_h_2, label='non mp')
    plt.legend()
    plt.grid(True)
    plt.xlabel('$n')
    plt.ylabel('$E[n]$') 

        
 

    plt.show()
    
    
    
    
########################
# numerically deriving the frequency response out of poles and zeros
# motivated by graphical relation of frequency response and pole-zero diagram
########################    
def my_freq_response(zer, pol, freqs, H_0 = 1):
    """
    gives frequency response of a discrete LTI system when zeros, poles and gain are provided
    
    IN: zeros, poles, H_0, number of freq. points or vector of freq samples
    
    OUT: H(Omega) at given freqs
    """
    # freq. vector
    if len(freqs)==1:
        num_freq_points = freqs
        Omega = np.linspace(-np.pi, np.pi, num_freq_points)
    else:
        num_freq_points = len(freqs)
        Omega = freqs
    
    # initialize freq. response
    H = np.zeros(num_freq_points, dtype = complex)
    
    # determine freq. response by distance and angle to zeros and poles
    for k in np.arange(num_freq_points):
        H[k] = np.prod( np.exp(1j*Omega[k])- zer) / np.prod( np.exp(1j*Omega[k])- pol )
    
    return freqs, H_0 * H


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
    
    #plt.plot( p.real, p.imag, '^', markersize=10)
    
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
    phi = np.linspace(0, 2*np.pi, 256)
    unit_circle = np.exp(1j*phi)
    for u in unit_circle:
        u_r.append(u.real)
        u_i.append(u.imag)
    plt.plot(u_r,u_i)



########################
# make it executable
########################
if __name__ == "__main__":
    main()
