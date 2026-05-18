# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 08:14:35 2014

@author: jaekel
"""


########################
# determine pswf functions by solving the corresponding eigen-equation
########################

# importing relevant stuff
import numpy as np
from scipy import signal

import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    # parameters: varying time-bandwidth-product
    # thereby fixing sample time to equal 1
    T = 1.0    
    t = np.linspace(-T, T, 1000)   
    dt = t[2]-t[1]
    
    TW_product = np.linspace(.01, 5.01, 10)
    W = TW_product / T / 4.0

    # initialize     
    kernel_matrix = np.zeros([len(t), len(t), len(W)], dtype=float)        
    Lambda_matrix = np.zeros([len(t), len(W)], dtype=complex)        
    X_matrix = np.zeros([len(t), len(t), len(W)], dtype=complex)        
    
    for w in np.arange(0,len(W)):
        for q in np.arange(0, len(t)):
                kernel_matrix[q,:,w] = own_sinc( t-t[q], W[w])

        Lambda_matrix[:, w], X_matrix[:,:,w] = np.linalg.eig(kernel_matrix[:,:,w])
        
#        # sorting 
        indices = (Lambda_matrix[:,w]).argsort()[::-1]
        Lambda_matrix[:,w] = Lambda_matrix[indices,w]
        X_matrix[:, :, w] = X_matrix[indices, :, w]
        
        

    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    
    plt.figure(14)    
    for w in np.arange(0,len(W)):
        plt.plot(np.arange(0,len(t)), dt*Lambda_matrix[:,w], label=str(W[w]))
    
    plt.ylabel('$\lambda_i$')
    plt.xlabel('$i$')    
    plt.grid(True)
    plt.title('Sorted Eigenvalues')
    #plt.legend()
    
    plt.figure(15)    
    #for w in np.arange(0,3):
    #    plt.plot(t, X_matrix[:,w,1], label='$\Psi(t)$')
    plt.plot(t, X_matrix[0,:,0])
    plt.ylabel('$\Psi_i(t)$')
    plt.xlabel('$t$')    
    plt.grid(True)
    plt.title('PSWF')
    #plt.legend()
    
    plt.show()    
    
    plt.show()
 
 
######################## 
# simplified sinc function avoiding NaN
########################
def own_sinc(x, w):
    x = np.asanyarray(x)
    y = np.pi*np.where(x==0, 1.0e-20, x)
    return np.sin(w*y) / y


    
########################
# make it executable
########################
if __name__ == "__main__":
    main()