# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:00:08 2014

@author: jaekel
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math



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
# main function
########################
def main():
    
    # set time resp. pulse interval and related parameters
    t_min=-0.0
    t_max=10.0
    t_s= 0.01                                 # sample time
    t=np.arange(t_min, t_max+t_s, t_s)
    
    f0=1 
    t_sample=.1
    t_sampled=np.arange(t_min, t_max+t_sample, t_sample)
    
    # frequency regime
    #f=np.arange( -1/(2*t_s), 1/(2*t_s)+1/(t_max-t_min), 1/(t_max-t_min))
    
    # original signal and sampled version thereof 
    x=np.sin(2*np.pi*f0*t)
    
    x_sampled=np.sin(2*np.pi*f0*t_sampled)
    
    # for illustration quantization of the signal may be included
    # w: number of bits corresponding to 2^w levels
    # input is assumed to be within [-1, 1]
    # define parameters
    w = 3
    quant = my_quantizer( w )

    x_quant = quant.do_quantization(x_sampled )



    x_recon=0*t
    x_recon_quant=0*t
    
    for k in np.arange(0, len(t_sampled)):

        argument = np.pi/t_sample*(t-k*t_sample)
        recon_sum_part = np.sin( argument ) / argument
        positions_nan = np.isnan(recon_sum_part)
        recon_sum_part[ positions_nan ] = 1
        
        x_recon += x_sampled[k] * recon_sum_part
        x_recon_quant += x_quant[k] * recon_sum_part        


    # plotting
    font = {'size'   : 36}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
       

    plt.figure(1)
    
    plt.subplot(121)
    plt.plot(t, x, label='$x(t)$')
    plt.plot(t_sampled, x_sampled,'ro', label='$x[n]$')
    plt.plot(t_sampled, x_quant,'g^', markersize=10, label='$x_q[n]$')
    plt.grid(True); plt.xlabel('$t/\mathrm{s}$');  #plt.ylabel('$x(t)$')
    plt.legend(loc='upper right')
    plt.axis(xmin=0, xmax=1.5, ymin=-1.1, ymax=1.1)
    plt.title('Signal und Samples')

    plt.subplot(122)
    plt.plot(t,x_recon, label='$x(t)$')
    plt.plot(t, x_recon_quant, label='$x_q(t)$')
    plt.grid(True); plt.xlabel('$t/\mathrm{s}$');  #plt.ylabel('$x(t)$')
    plt.axis(xmin=0, xmax=1.55, ymin=-1.1, ymax=1.1)
    plt.title('Rekonstruiertes Signal')
    plt.legend(loc='upper right')    
    
    plt.show()
    
  
  
########################
# make it executable
########################
if __name__ == "__main__":
    main()