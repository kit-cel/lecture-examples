# -*- coding: utf-8 -*-

########################
######## importing
########################
  
from scipy import signal

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


########################
# quantizer class
########################
class my_quantizer():
    '''
    Class for quantizing signals
    
    Parameters:     w: bits per quant. level --> 2^w levels
                    
    NOTE: quantization is done as midtread and using symmetrical amplitudes,
            thus discardong one quant. level/word

    NOTE: quantization is assuming that input values are within [-x_max, x_max]
            default is x_max = 1

    Methods:        __init__
                    find_levels: returns quant. levels based on midtread approach
                        decision_bounds = ((levels + np.roll(levels, -1))/2)[:numb_levels-1]

    '''

        
    #####################            
    def __init__( self, w, x_max = 1):
        ''' 
        initialze object
        '''

        self.w = w
        self.K = 2**w
        self.x_max = x_max
        self.Delta = 2 * x_max / self.K

        self.find_levels()


    #####################               
    def find_levels( self ):
        '''
        find quantization levels and symmetrical decision boundaries
        '''

        levels = [ 0. ]
        for k in range( 1, 2**( self.w - 1 ) ):
            levels.extend( [ -k * self.Delta, k * self.Delta ] )
            
        self.levels =  np.array( np.sort( levels ) )
        self.decision_bounds = ( ( self.levels + np.roll( self.levels, - 1 ) ) / 2 ) [ : self.K - 2 ]


    #####################               
    def do_quantization( self, x ):
        '''
        find quantization levels and assign them to signal x_q
        '''

        x_q = []
        for xi in x:

            bigger = np.array( np.where( xi >= self.decision_bounds ) )
            
            if bigger.size:
                x_q.append( self.levels[np.max(bigger) +1 ] )
            else:
                x_q.append( self.levels[0])
        
        return np.array( x_q )


########################
# sample signal
########################
# get random function by randomly cascading 4 functions 
def get_function( x ):

    # get sum of two sine functions
    N_max = 15
    N = 3 + np.random.randint( N_max )

    coeffs = np.random.rand( N ) * ( -1 )**np.random.randint( 2, size=N )
    poly = np.polynomial.polynomial.Polynomial( coeffs )

    functions = [ np.sin, np.arcsin, np.sqrt, np.log, np.exp, np.i0, np.sinh, np.arcsinh, np.sinc, poly]

    fct = np.random.choice( functions, 4)
    y = x
    for f in fct:
        y = f( y )
        y[ np.isnan( y ) ] = 0
        y[ np.isinf( y ) ] = 0

    return y


########################
# main function
########################
def main():
    
    # define parameters
    w = 3
    quant = my_quantizer( w ) 


    # dft length
    N = 32
    t_s  = 0.125
    t_array = np.arange( 0, N ) * t_s

    n_array = np.arange(0,N)
    k_array = np.arange(0,N)

    # multiples of signal length after zero padding
    M = 100
    k_padded=np.arange(0, N, 1.0/M)

    # parameter from the slides     
    ell = 5.0
    delta_ell=.25

    f_0 = ( ell + delta_ell ) / N / t_s 
    
    # signals (being harmonics)
    x = np.exp( 1j * 2 * np.pi * f_0 * t_array ) 
    x = get_function( t_array ) 

    x_max = np.max( np.abs( x ) )
    X = np.fft.fft( x )

    x_padded = np.append( x, np.zeros((M-1)*N))
    X_padded = np.fft.fft( x_padded )
    
    # find quantized versions
    x_quant = quant.do_quantization( x / x_max ) * x_max
    X_quant = np.fft.fft( x_quant ) 

    x_quant_padded = np.append( x_quant, np.zeros((M-1)*N))
    X_quant_padded = np.fft.fft( x_quant_padded )  


    # plotting
    font = {'size'   : 36}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    
    plt.close('all')

    plt.figure(1, figsize=(16, 9) )
    plt.subplot(221)
    plt.stem(n_array, x.real,label='$Re\{ x[n]\}$')
    plt.grid(True);  #plt.xlabel('$n$'); #plt.ylabel('$Re\{ x[n]\}$')
    plt.legend(loc='upper right')
    
    plt.subplot(222)        
    plt.stem(k_array, np.abs(X), label='$|X[k]|$')
    plt.plot(k_padded, np.abs( X_padded ),'--')
    plt.grid(True);  #plt.xlabel('$k$'); #plt.ylabel('$|X[k]|$')
    plt.legend(loc='upper right')    
    
    plt.subplot(223)
    plt.stem(n_array, x_quant.real,label='$Re\{ x_q[n]\}$')
    plt.grid(True);  plt.xlabel('$n$'); #plt.ylabel('$Re\{ x[n]\}$')
    plt.legend(loc='upper right')
    
    plt.subplot(224)        
    plt.stem(k_array, np.abs(X_quant), label='$|X[k]|$')
    plt.plot(k_padded, np.abs( X_quant_padded ),'--')
    plt.grid(True);  plt.xlabel('$k$'); #plt.ylabel('$|X[k]|$')
    plt.legend(loc='upper right')    

    plt.show()
    

   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
