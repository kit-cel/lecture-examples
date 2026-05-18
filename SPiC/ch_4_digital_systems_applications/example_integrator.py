# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:10:50 2015

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
    
    K1 = 10
    K2 = 30

    Omega = np.arange(0, np.pi, np.pi/1024)        
    
    h1 = np.ones(K1)/np.sqrt(K1)
    h2 = np.ones(K2)/np.sqrt(K2)
    
    H1_num = np.fft.fft( np.append( h1, np.zeros(1024-K1) ) )
    H2_num = np.fft.fft( np.append( h2, np.zeros(1024-K2) ) )  
    
    H1_an = np.exp(-1j*Omega/2*(K1-1)) * np.sin(Omega*K1)/np.sin(Omega)
    H1_an[ np.isnan(H1_an)] = K1
    H1_an *= np.max(H1_num)/np.max(H1_an)


    # plotting
    plt.figure()
    plt.plot(Omega, np.abs(H1_num), label='$num., K=$'+str(K1))
    plt.plot(Omega, np.abs(H1_an), 'o', label='$an., K=$'+str(K1))    
    
    plt.plot(Omega, np.abs(H2_num), label='$num., K=$'+str(K2))    
    plt.legend(loc = 'upper right')
    plt.grid(True)
    plt.xlabel('$\\Omega$')
    plt.ylabel('$|H(\\Omega)|$')

    plt.show()
    
    # that's it folks
    print('Done!')
   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
