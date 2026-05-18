# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 11:25:31 2015

@author: jaekel
"""

from __future__ import division
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

   
import platform
import sys

if platform.system() == "Windows":
    sys.path.append('Z:\\Simulationen_Programme/Kommunikationssysteme')
    
if platform.system() == "Linux":   
    sys.path.append('/home/jaekel/Dokumente/Simulationen_Programme/Kommunikationssysteme/')
        
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
    
    # filtering with a windowed rrc
    rrc_win = rrc*uc.windowing(rrc, 'Gauss')
    s_filt_rrc_win = np.convolve(s, rrc_win)
    
    # at Rx
    s_mf = np.convolve( s_filt_rrc, rrc)
    s_mf_win = np.convolve(s_filt_rrc_win, rrc_win)

    # frequency regime
    S_mf = np.fft.fft( np.append( s_filt_rrc, np.zeros(9*len(s_filt_rrc  ) ) ) )
    S_mf_win = np.fft.fft( np.append( s_filt_rrc_win, np.zeros(9*len(s_filt_rrc_win  ) ) ) )

    Omega = np.arange( 0, 2*np.pi, 2*np.pi/len(S_mf))
    
    # plotting
    plt.figure()
    plt.subplot(121)
    plt.plot( np.arange(np.size(s_filt_rrc)), np.real(s_filt_rrc),'o', label='RRC')
    plt.plot( np.arange(np.size(s_filt_rrc_win)), np.real(s_filt_rrc_win),'x', label='windowed RRC')    
    plt.grid(True); plt.xlabel('$n$');  plt.title('after Tx RRC'); plt.legend(loc='upper right')

    plt.subplot(122)
    plt.plot( np.arange(np.size(s_mf)), np.real(s_mf),'o', label='RRC')
    plt.plot( np.arange(np.size(s_mf_win)), np.real(s_mf_win),'x', label='windowed RRC')    
    plt.grid(True); plt.xlabel('$n$');  plt.title('after Rx MF'); plt.legend(loc='upper right')
   



    plt.figure()
    plt.subplot(121)
    plt.plot( Omega, np.abs(S_mf), label='$|$RRC$|$')
    plt.plot( Omega, np.abs(S_mf_win), label='$|$windowed RRC$|$' )    
    plt.grid(True)
    plt.xlabel('$\\Omega$')
    plt.legend(loc='upper right')
    
    plt.subplot(122)
    plt.plot( Omega, np.angle(S_mf), label='$\\angle$ RRC')
    plt.plot( Omega, np.angle(S_mf_win), label='$\\angle$ windowed RRC' )    
    plt.grid(True)
    plt.xlabel('$\\Omega$')    
    plt.legend(loc='upper right')
     
    plt.show()
    
    # that's it
    print('Done!')
 
 
 

        
       
       
########################
# make it executable
########################
if __name__ == "__main__":
    main()