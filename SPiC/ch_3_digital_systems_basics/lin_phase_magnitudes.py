# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:10:55 2015

@author: jaekel
"""

from __future__ import division
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

   
########################
# parameters
########################

# plotting options to have "good looking" figures
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
       

  
########################
# main function
########################
def main():
    
    h_u_s = [.1, .15, .2, .3, 1, .3, .2, .15, .1]
    h_u_s /= np.linalg.norm(h_u_s)
    
    h_g_s = [.1, .15, .3, .7, .7, .3, .15, .1]    
    h_g_s /= np.linalg.norm(h_g_s)
    
    h_u_a = [.1, .2, .4, .6, 0, -.6, -.4, -.2, -.1]    
    h_u_a /= np.linalg.norm(h_u_a)
    
    h_g_a = [.1, .2, .4, .7, -.7, -.4, -.2, -.1]    
    h_g_a /= np.linalg.norm(h_g_a)
  
    N_fft = 1024
    Omega = np.arange( -np.pi, np.pi, 2*np.pi/N_fft)
  
    H_u_s = get_H(h_u_s, N_fft)
    H_u_a = get_H(h_u_a, N_fft)    
    H_g_s = get_H(h_g_s, N_fft)    
    H_g_a = get_H(h_g_a, N_fft)        
    
   
    # plotting
    plt.figure()
    
    # impulse responses
    plt.subplot(4,2,1)
    plt.stem( np.arange(np.size(h_u_s)), h_u_s,'ob-')
    plt.grid(True); plt.title('$h[k]$'); plt.ylabel('uL-s')
    plt.axis([0, 10, -1, 1])
    
    plt.subplot(4,2,3)
    plt.stem( np.arange(np.size(h_g_s)), h_g_s,'ob-')
    plt.grid(True); plt.ylabel('gL-s')
    plt.axis([0, 10, -1, 1])    
    
    plt.subplot(4,2,5)
    plt.stem( np.arange(np.size(h_u_a)), h_u_a,'ob-')
    plt.grid(True); plt.ylabel('uL-a')
    plt.axis([0, 10, -1, 1])
        
    plt.subplot(4,2,7)
    plt.stem( np.arange(np.size(h_g_a)), h_g_a,'ob-')
    plt.grid(True); plt.xlabel('$k$'); plt.ylabel('gL-a')
    plt.axis([0, 10, -1, 1])
    
    # magnitudes
    plt.subplot(4,2,2)
    plt.plot( Omega, np.abs(H_u_s))
    plt.grid(True); plt.title('$|H(\Omega)|$')
    plt.axis([-np.pi, np.pi, 0, 2.5])

    plt.subplot(4,2,4)
    plt.plot( Omega, np.abs(H_g_s))
    plt.grid(True); 
    plt.axis([-np.pi, np.pi, 0, 2.5])
    
    plt.subplot(4,2,6)
    plt.plot( Omega, np.abs(H_u_a))
    plt.grid(True); 
    plt.axis([-np.pi, np.pi, 0, 2.5])
    
    plt.subplot(4,2,8)
    plt.plot( Omega, np.abs(H_g_a))
    plt.grid(True); plt.xlabel('$\Omega$');         
    plt.axis([-np.pi, np.pi, 0, 2.5])
        
              
    
    plt.show()
    
    # that's it
    print('Done!')
 
 
########################
# for easy of notation
######################## 
def get_H(h, N_fft):
    return np.fft.fftshift( np.fft.fft( h, N_fft ) ) 
         
       
       
########################
# make it executable
########################
if __name__ == "__main__":
    main()