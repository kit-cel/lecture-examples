# -*- coding: utf-8 -*-
"""
Created on Wed Mar 04 15:41:59 2015

@author: jaekel
"""

from __future__ import division
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

   
import platform
import sys

if platform.system() == "Windows":
    sys.path.append('Z:\\Simulationen_Programme/kommunikationssysteme')
    
if platform.system() == "Linux":   
    sys.path.append('/home/jaekel/Dokumente/Simulationen_Programme/kommunikationssysteme/')
        
import utils_communications as uc   



 
########################
# parameters
########################

n_symb = 1

# set time resp. pulse interval and related parameters
t_min = 0.0
t_max = float(n_symb)
t_s = 0.1                                 # sample time
t = np.arange(t_min, t_max+t_s, t_s)

t_symb = 1.0 

# parameters for rrc
beta = 0.33

# plotting options to have "good looking" figures
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
       

  
########################
# main function
########################
def main():
    
    #####    
    # with "real rrc" 
    #####        
    n_up = 8            # samples per symbol
    syms_per_filt = 2  # symbols per filter (plus minus in both directions)
    
    K_filt = 2*syms_per_filt*n_up+1         # length of the fir filter

    print(K_filt) 
    
    #########################################        
    # for debugging: simple pulse        
    s = 1#np.append(1, np.zeros(n_symb-1))
    #########################################
        
        
    #####
    # apply rrc        
    #####
    rrc = uc.get_rrc_ir(n_up*syms_per_filt*2+1, n_up, t_symb, beta)
    E_rrc = np.sum(np.abs(rrc)**2)
    
    rrc = rrc/np.sqrt(E_rrc)

    # filtering without using a window    
    s_filt_rrc = np.convolve(s, rrc)
    
    S = np.fft.fftshift( np.fft.fft( np.append(s_filt_rrc, np.zeros(9*np.size(s_filt_rrc)))))
    S_mag = np.abs(S)
    S_angle = np.angle(S)
    S_group = -np.diff(S_angle)


    print( S_group )

    # filtering with a windowed rrc
    rrc_win = rrc*uc.windowing(rrc, 'Gauss')
    s_filt_rrc_win = np.convolve(s, rrc_win)
    
    S_win = np.fft.fftshift( np.fft.fft( np.append(s_filt_rrc_win, np.zeros(9*np.size(s_filt_rrc_win)))))
    S_mag_win = np.abs(S_win)
    S_angle_win = np.angle(S_win)
    S_group_win = -np.diff(S_angle_win)
    
    
    # plotting
    plt.figure()
    plt.subplot(131)
    plt.plot( np.arange(np.size(s_filt_rrc)), np.real(s_filt_rrc),'o', label='RRC')
    #plt.plot( np.arange(np.size(s_filt_rrc_win)), np.real(s_filt_rrc_win),'x', label='windowed RRC')    
    plt.grid(True); plt.xlabel('$n$');  plt.title('$RRC[n]$'); 
    #plt.legend(loc='upper right')
        
    plt.subplot(132)        
    f = np.linspace(-1/(2*t_s*n_up), 1/(2*t_s*n_up), num = np.size(S_mag))
    
    plt.plot(f, 10*np.log10(S_mag/np.max(S_mag)), label='RRC')
    #plt.plot(f, 10*np.log10(S_mag/np.max(S_mag_win)), label='windowed RRC')    
    plt.axis( [np.min(f), np.max(f), -50, 0])    
    plt.grid(True); plt.xlabel('$f/\mathrm{Hz}$'); plt.title('$|H(f)|$'); 
    #plt.legend( loc='lower left')
            
    plt.subplot(133)        
    plt.plot(f[:-1], np.remainder(S_group, 2*np.pi) - np.pi, label='RRC')
    #plt.plot(f[:-1], np.remainder(S_group_win, 2*np.pi) - np.pi, label='windowed RRC')    
    plt.axis( [np.min(f), np.max(f), -np.pi, np.pi])
    plt.grid(True); plt.xlabel('$f/\mathrm{Hz}$'); plt.title('$\\tau_\\mathrm{g}(f)$')
    #plt.legend(loc='upper right')
          
    plt.show()
    
    # that's it
    print('Done!')
 
 
 

        
       
       
########################
# make it executable
########################
if __name__ == "__main__":
    main()