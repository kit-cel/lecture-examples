# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 08:34:10 2014

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
    N = 21
    n = np.arange(N)
    
    a = .8

    #num = [1, 0]
    #denom = [1, -a]
    #[zeros, poles, H_0] = signal.tf2zpk( num, denom)
    
    zeros = np.array([0, 0, 0])
    poles = np.array([.5, .25+.8j, .25-.8j])
    num, denom = signal.zpk2tf(zeros, poles, 1)
    
    n_o, h = signal.dimpulse( [num, denom, 1])
    h = np.squeeze(np.array(h))
    
    freqs = np.linspace(-np.pi, np.pi, 512)    
    w, H = signal.freqz(num, denom, worN=freqs, whole=True)    
    H = np.squeeze(np.array(H))    
    

    # illustrate inverse system
    zeros_i = poles
    poles_i = zeros
    num_i, denom_i = signal.zpk2tf(zeros_i, poles_i, 1)
    
    n_o_i, h_i = signal.dimpulse( [num_i, denom_i, 1])
    h_i = np.squeeze(np.array(h_i))
    
    w_i, H_i = signal.freqz(num_i, denom_i, worN=freqs, whole=True)    
    H_i = np.squeeze(np.array(H_i))  
    
    
    # plotting
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)  

    plt.figure(1)
    plt.subplot(121)
    plot_zeros_poles(zeros, poles, [-1.1, 1.1, -1.1, 1.1], 'Re$\{z\}$', 'Im$\{z\}$')
    plt.annotate('z plane',
            xy=(0,.45),
            xytext=(-1,1),
            )  
    plt.annotate('unit circle',
            xy=(.3,.95),
            xytext=(.5,1),
            arrowprops={'facecolor':'green','shrink':0.05},            
            )  
    plt.annotate('$(3)$',
            xy=(0,.45),
            xytext=(0,0.1),
            )     
            
            
    plt.subplot(122)
    plt.stem(n_o, h.real, label='$h[n]$')#, n_o, h.imag, 'o')
    #plt.axis([-.1, 30.1, -.1, 1.1])
    plt.grid(True)        
    plt.xlabel('$n$')
    #plt.ylabel('$h[n]$')
    plt.legend(loc='upper right')
    
    
    plt.figure(2)    
    plt.subplot(121)
    plt.plot(freqs, abs(H), label='$|H(\\Omega)|$')
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    #plt.ylabel('$|H(\Omega)|$')
    plt.legend(loc='upper right')    

    plt.subplot(122)    
    plt.plot(freqs, np.angle(H), label='$\\angle H(\\Omega)$')
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    #plt.ylabel('$\phi(\Omega)$')    
    plt.legend(loc='upper right')


    plt.figure(3)
    plt.subplot(121)
    plot_zeros_poles(zeros_i, poles_i, [-1.1, 1.1, -1.1, 1.1], 'Re$\{z\}$', 'Im$\{z\}$')
    plt.annotate('z plane',
            xy=(0,.45),
            xytext=(-1,1),
            )  
    plt.annotate('unit circle',
            xy=(.3,.95),
            xytext=(.5,1),
            arrowprops={'facecolor':'green','shrink':0.05},            
            )  
    plt.annotate('$(3)$',
            xy=(0,.45),
            xytext=(0,0.1),
            )     
            
            
    plt.subplot(122)
    plt.stem(n_o, h_i.real,label='$h[n]$')#, n_o, h.imag, 'o')
    #plt.axis([-.1, 30.1, -.1, 1.1])
    plt.grid(True)        
    plt.xlabel('$n$')
    #plt.ylabel('$h[n]$')
    plt.legend(loc='upper right')


    plt.figure(4)
    plt.subplot(121)
    plt.plot(freqs, abs(H_i), label='$|H(\\Omega)|$')
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    #plt.ylabel('$|H_i(\Omega)|$')
    plt.legend(loc='upper right')

    plt.subplot(122)    
    plt.plot(freqs, np.angle(H_i), label='$\\angle H(\\Omega)$')
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    #plt.ylabel('$\phi_i(\Omega)$')  
    plt.legend(loc='upper right')

    
    plt.figure(5)
    plt.subplot(121)
    plt.plot(freqs, abs(H)*abs(H_i))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$|H(\Omega)|\cdot |H_i(\Omega)|$')

    plt.subplot(122)    
    plt.plot(freqs, np.angle(H)+np.angle(H_i))
    plt.grid(True)
    #plt.axis([-np.pi-.1, np.pi+.1, -.1, 5.1])
    plt.xlabel('$\Omega$')
    plt.ylabel('$\phi(\Omega)+\phi_i(\Omega)$')  
    
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
