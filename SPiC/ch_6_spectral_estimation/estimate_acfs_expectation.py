# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 09:05:48 2014

@author: jaekel
"""

########################
# illustrating asymptotically unbiased estimation of psd
########################

# importing relevant stuff
import numpy as np
import own_utils_est_psd as own_spec

import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################
def main():
    # parameters: number of points in first resp. second function and in frequency domain
    N_1 = int( 1e3 )
       
    N_acf_1 = np.arange(-N_1+1, N_1, 1)
    
    # number of freq. points and freq. range
    N_freq = 512            
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
    # number of realizations for averaging    
    N_aver = int( 1e2 )
    
    ###
    # first function
    ###
    # relative width of rectangular 
    relative_width_rect = 0.3
    N_rect = int(N_1*relative_width_rect)
    f = np.append( np.ones( N_rect), np.zeros( N_1-N_rect))
    
    
    # loop for realizations
    acf_biased_1 = np.zeros(len(N_acf_1))
    acf_unbiased_1 = np.zeros(len(N_acf_1))
    
    psd_p_1 = np.zeros(len(Ome))
    psd_c_biased_1 = np.zeros(len(Ome))
    psd_c_unbiased_1 = np.zeros(len(Ome))
    
    n=0
    while n < N_aver:
    
        # first function    
        #f_1 = f + np.random.normal(0.0, 0.2, N_1)
        f_1 = np.random.normal(0.0, 1.0, N_1)
        #f_1 = np.sin(1.5*np.arange(0, N_1)) + np.random.normal(0.0, 0.2, N_1)
        
        acf_biased_1 = 1./(n+1) *(float(n)*acf_biased_1 + own_spec.est_acf(f_1, 'biased') )
        acf_unbiased_1 = 1./(n+1) *(float(n)*acf_unbiased_1 + own_spec.est_acf(f_1, 'unbiased') )
        
     
        # find periodogram by simple fft and abs()**2
        psd_p_1 = 1./(n+1) *(n*psd_p_1 + own_spec.find_periodogram(f_1, Ome) )
        
        psd_c_biased_1 = 1./(n+1) *(n*psd_c_biased_1 + own_spec.find_correlogram( acf_biased_1, Ome) )
        psd_c_unbiased_1 = 1./(n+1) *(n*psd_c_unbiased_1 + own_spec.find_correlogram( acf_unbiased_1, Ome) )
   
        n += 1
    
        # show progress
        done = float(n)/N_aver*100.0
        print('Done: %3.2f percent' % done)
     
     
     
     
    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    
    # plot acfs
    plt.figure(1)    
    plt.subplot(111)
    plt.plot(N_acf_1, acf_unbiased_1,'b', label='$L=N-k$')      
    plt.plot(N_acf_1, acf_biased_1,'r', label='$L=N$')
    #plt.plot(N_acf_1, (1-abs(N_acf_1)/N_1)*acf_wgn_unbiased_1,'--')
    plt.xlabel('$k$')
    plt.ylabel('$\hat{r}[k]$')   
    plt.grid(True)    
    plt.title('$N=$'+str(N_1))
    plt.legend(loc='upper right') 
   
    
    # plot psd1
    plt.figure(2)    
    
    plt.subplot(131)    
    plt.plot(Ome, psd_p_1)      
    #plt.plot(Ome_corr_1, psd_biased_1, label='corr., b')          
    #plt.title('$N=$'+str(N_1))    
    plt.grid(True); #plt.legend(loc='upper right') 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_p(\Omega)$')   
    plt.axis([-4, 4, 0, 2])
    
    plt.subplot(132)    
    plt.plot(Ome, psd_c_biased_1, label='$L=N$')      
    plt.grid(True); plt.legend(loc='upper right') 
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')       
    plt.axis([-4, 4, 0, 2])
    
    plt.subplot(133)        
    plt.plot(Ome, psd_c_unbiased_1, label='$L=N-k$')             
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.grid(True); plt.legend(loc='upper right') 
    plt.axis([-4, 4, 0, 2])
    
    
    plt.show()
        



    
########################
# make it executable
########################
if __name__ == "__main__":
    main()