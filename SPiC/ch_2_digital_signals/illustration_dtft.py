# -*- coding: utf-8 -*-
"""
Created on Wed Apr 01 16:30:12 2015

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
    t_min=-10.0
    t_max=10.0    
      
    # sampled version 
    t_s = .25
    
    t_samp = np.arange(t_min, t_max, t_s)
    x_samp = (signal.gausspulse(t_samp, 1, retenv=1))[1] 
    x_samp = x_samp * np.sqrt( 1 / sum(x_samp**2) )

    # determine according Nyquist band and fft
    f_samp = np.arange( -1/(t_s*2.0), 1/(t_s*2.0), 1/(t_max-t_min) )    
    X_samp = np.fft.fftshift(np.fft.fft(x_samp))
    
    # upsampling
    M = 4
    t_s_up = t_s/M
    t_up = np.arange(t_min, t_max, t_s_up)
    
    # upsample signal by factor M using lambda function
    x_up = np.zeros( M * len( x_samp ) )
    x_up[ :: M ] = x_samp
    
    x_cont = (signal.gausspulse(t_up, 1, retenv=1))[1] 
    x_cont /= (np.max(x_cont)/np.max(x_samp))

    # determine according Nyquist band and fft
    f_up = np.arange( -1/(t_s_up*2.0), 1/(t_s_up*2.0), 1/(t_max-t_min) )    
    X_up = np.fft.fftshift(np.fft.fft(x_up))
    
  
    
    
    # plotting
    font = {'size'   : 36}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    
   
    
    plt.figure(1)
    plt.subplot(121)
    plt.plot(t_samp, x_samp, 'o', markersize=10, label='$x[n]$')
    plt.plot(t_up, x_cont, '', label='$x(t)$')
    plt.grid(True)
    plt.axis(xmin=-3, xmax=3, ymin=-.1, ymax=.5)    
    plt.xlabel('$t/\mathrm{s}$')
    #plt.ylabel('$x(t), x[n]$')
#    plt.annotate('Samples',
#            xy=(0,.435),
#            xytext=(1,0.45),
#            arrowprops={'facecolor':'green','shrink':0.05},
#            )    
    plt.legend( loc='upper right' )
    
    plt.subplot(122)        
    plt.plot(f_up, abs(X_up)**2, '', label='$|\mathcal{F}\{x[n]\}(f)|^2$')
    plt.plot(f_samp, np.abs(X_samp)**2, '', label='$|X(f)|^2$')

    plt.grid(True)
    plt.xlabel('$f/\mathrm{Hz}$')
    #plt.ylabel('$|X(f)|^2$')
    plt.annotate('Nyquist Spektrum',
            xy=(0,10.5),
            xytext=(1.5,11),
            arrowprops={'facecolor':'green','shrink':0.05},
            )    
    plt.annotate('Period. Wdhl.',
            xy=(-4,10.5),
            xytext=(-7,11),
            arrowprops={'facecolor':'green','shrink':0.05},
            )      
    plt.legend( loc='upper right' )
    
    plt.show()
    
    # that's it folks
    print('Done!')
   
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
