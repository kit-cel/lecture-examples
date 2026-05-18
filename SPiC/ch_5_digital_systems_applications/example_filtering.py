# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 13:36:56 2015

@author: jaekel
"""


########################
######## importing and stuff
########################
  
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)


########################
# main function
########################
def main():


    # parameters
    length_x = int( 1e4 )

    K = 20    
    length_h = 50
    
    length_fft = 512
    
    length_p = 100
    
    # generate white noise as input signal
    t_s = 1.0
    t = np.arange(length_x)

    x = np.ones(len(t)) + np.sin(2*np.pi*1.0/4*t) + np.sqrt(1)*np.random.randn(length_x)
    
    # define several impulse responses
    h_1 = np.zeros(length_h)
    h_1[:K] = 1
    h_1 /= np.linalg.norm(h_1)
    
    h_2 = np.zeros(length_h)    
    h_2[:K] = np.sign( np.sin(2*np.pi*1.0/6*np.arange(K))) 
    h_2 /= np.linalg.norm(h_2)
    
    h_3 = np.zeros(length_h)
    h_3[:K] = (-1)**np.arange(K)
    h_3 /= np.linalg.norm(h_3)
   
   
    # get outputs
    y_1 = np.convolve(x, h_1)
    y_2 = np.convolve(x, h_2)    
    y_3 = np.convolve(x, h_3)
    
    # spectra
    f = np.linspace(-1/(2*t_s), 1/(2*t_s), length_fft)
    
    X = find_periodogram(x, 2*np.pi*f)
    Y_1 = find_periodogram(y_1, 2*np.pi*f) 
    Y_2 = find_periodogram(y_2, 2*np.pi*f) 
    Y_3 = find_periodogram(y_3, 2*np.pi*f)     
    
    
    # plotting
    x_plot = x[:length_p]
    y_1_plot = y_1[:length_p]
    y_2_plot = y_2[:length_p]    
    y_3_plot = y_3[:length_p]    
    
    plt.figure(1)
    plt.clf()

    plt.subplot(331)
    plt.plot( np.arange(len(h_1)), h_1, 'o')
    plt.grid(True)
    #plt.axis([0, length_h, -1.1, 1.1])
    plt.title('$h[n]$')

    plt.subplot(332)

    plt.plot( np.arange(len(x_plot)), x_plot, label='$x[n]$')    
    plt.plot( np.arange(len(y_1_plot)), y_1_plot, 'r', label='$y[n ]$')
    plt.grid(True)
    #plt.legend(loc='upper right')
    plt.title('$x[n], y[n]$')
    #plt.axis([0, length_p, -1.1, 10])   
    plt.xlim( xmax=length_p)    

    plt.subplot(333)
    plt.plot(f, X, label='$|X(f)|^2$')     
    plt.plot(f, Y_1, 'r', label='$|Y(f)|^2$') 
    plt.grid(True)
    #plt.legend(loc='upper right')        
    plt.title('$|X(f)|^2, |Y(f)|^2$')

    plt.subplot(334)
    plt.plot( np.arange(len(h_2)), h_2, 'o')
    plt.grid(True)
    #plt.axis([0, length_h, -1.1, 1.1])
    plt.xlabel('$n$')     
    
    plt.subplot(335)
    plt.plot( np.arange(len(x_plot)), x_plot, label='$x[n]$')    
    plt.plot( np.arange(len(y_2_plot)), y_2_plot, 'r', label='$y[n ]$')
    plt.grid(True)
    #plt.legend(loc='upper right')
    #plt.axis([0, length_p, -1.1, 10])
    plt.xlim( xmax=length_p)
    plt.xlabel('$n$') 
    
    plt.subplot(336)
    plt.plot(f, X, label='$|X(f)|^2$')     
    plt.plot(f, Y_2, 'r', label='$|Y(f)|^2$') 
    plt.grid(True)
    #plt.legend(loc='upper right')    
    plt.xlabel('$f/\mathrm{Hz}$')    

    plt.subplot(337)
    plt.plot( np.arange(len(h_3)), h_3, 'o')
    plt.grid(True)
    #plt.axis([0, length_h, -1.1, 1.1])
    plt.xlabel('$n$')     
    
    plt.subplot(338)
    plt.plot( np.arange(len(x_plot)), x_plot, label='$x[n]$')    
    plt.plot( np.arange(len(y_3_plot)), y_3_plot, 'r', label='$y[n ]$')
    plt.grid(True)
    #plt.legend(loc='upper right')
    #plt.axis([0, length_p, -1.1, 10])
    plt.xlim( xmax=length_p)
    plt.xlabel('$n$') 
    
    plt.subplot(339)
    plt.plot(f, X, label='$|X(f)|^2$')     
    plt.plot(f, Y_3, 'r', label='$|Y(f)|^2$') 
    plt.grid(True)
    #plt.legend(loc='upper right')    
    plt.xlabel('$f/\mathrm{Hz}$')  
    
     
   
    plt.show()
    
    # that's it folks
    print('Done!')
   
    
    
########################
# own periodogram estimator
########################
def find_periodogram(y, omega):
    """
    estimates periodogram out of the given observation at the frequencies specified in omega
    
    IN: observation y, frequencies
    OUT: psd estimator
    """
    N = len(y)
    per = np.zeros(len(omega), dtype=complex) 
        
    for p in np.arange(0, N):
        per += y[p]*np.exp(-1j*omega*(p+1))
        
    per = (abs(per)**2)/N
        
    return per  
    

########################
# make it executable
########################
if __name__ == "__main__":
    main()
