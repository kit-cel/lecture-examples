# -*- coding: utf-8 -*-
"""
Created on Wed May  7 11:08:55 2014

@author: jaekel
"""

########################
######## importing
########################
  
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt



########################
# main function
########################
def main():
    
    # dimension of the mimo system
    N_tx = 4    
    N_rx = 128
    
    # number of trials for estimation
    numb_trials = 1e3
    
    # list storing the ranks
    ranks = np.zeros(numb_trials)
    
    for n in np.arange(0, numb_trials):
        # find random channel matrix of flat block fading with circular gaussians
        H = 1./np.sqrt(2) * ( np.random.normal(size=(N_rx, N_tx)) + 1j*np.random.normal(size=(N_rx, N_tx)) )        

        # find rank and save
        ranks[n] = np.linalg.matrix_rank( H )
        
        print float(n)/numb_trials*100, '% done'

       
    # determine histogram
    hist, bins = np.histogram( ranks, bins = np.arange(1,max(N_tx, N_rx)+2), density = True)
    
    # plotting
    plt.figure(18)
    bar_width = .8
    plt.bar( bins[:-1]-bar_width/2, hist, width=bar_width)
    plt.xlim([min(bins)-1, max(bins)-1])
    plt.grid(True)
    plt.show()
    
    
########################
# make it executable
########################
if __name__ == "__main__":
    main()
    