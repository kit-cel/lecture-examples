# -*- coding: utf-8 -*-
"""
Created on Fri May 20 08:12:34 2016

@author: jaekel
"""

from __future__ import division

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
   
   
########################
# parameters
########################


# number of symbols per OFDM symbol   
N_ofdm = 8

# set time resp. pulse interval and related parameters
T_symb = 1.0 

T_ofdm = T_symb * N_ofdm
delta_f = 1 / T_ofdm


# plotting options to have "good looking" figures
font = {'size'   : 30}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
       
plt.close('all')
  
########################
# main function
########################
def main():
    
    # specify frequency-selective channel
    L = 4    
    h = [1, .5, .25, .25]

    # extend to have same length as other signals for FFT    
    h_app = np.append( h, np.zeros( N_ofdm - L ) )    
    H = np.fft.fft( h_app )
    
    # generate vector to be transmitted;
    # choosing increasing integers for simplicity
    d = range( N_ofdm ) 
    
    x = np.fft.ifft( d )
    
    # apply zero prefix (_zp) and cyclic prefix (_cp)    
    s_zp = np.append( np.zeros( L ), x )
    s_cp = np.append( x[ N_ofdm - L : ], x )
    
    # received signal by convolution in time domain
    r_zp = np.convolve( s_zp, h )
    r_cp = np.convolve( s_cp, h )

    y_zp = r_zp[ L : L + N_ofdm ]
    y_cp = r_cp[ L : L + N_ofdm ]

    # back to symbol domain
    Y_zp = np.fft.fft( y_zp )
    
    Y_cp = np.fft.fft( y_cp )

    np.set_printoptions(suppress=True)
    np.set_printoptions( precision = 2 )

    print('d = {}\n'.format( d ) )
    
    print('H = {}\n'.format( H ) )

    print('Y_ZP = {}\n'.format( Y_zp ) )    
    print('Y_CP = {}\n'.format( Y_cp ) )    
    
    

    # plotting    

    # that's it folks    
    print( 'Done!' ) 
 
 
 
    
    
########################
# make it executable
########################
if __name__ == "__main__":
    main()