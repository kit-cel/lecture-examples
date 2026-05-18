# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 09:26:31 2014

@author: jaekel
"""


########################
# illustrating estimation of acf and psd
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
    N_1 = int( 1e2 )
    N_2 = int( 1e2 )
    
    N_freq = 512

    N_acf_1 = np.arange( - N_1 + 1, N_1, 1 )
    N_acf_2 = np.arange( - N_2 + 1, N_2, 1 )    
    
    
    ###
    # first function
    ###
    # relative width of rectangular 
    relative_width_rect = 0.3
    N_rect = int(N_1*relative_width_rect)
    f_1 = np.append( np.ones( N_rect), np.zeros( N_1-N_rect))
    f_1 += np.random.normal(0.0, .2, N_1)
    #f_1 = np.exp(1j*1.5*np.arange(0, N_1))

    acf_biased_1 = own_spec.est_acf(f_1, 'biased')
    acf_biased_1 = acf_biased_1 / np.max(np.abs(acf_biased_1))
    
    acf_unbiased_1 = own_spec.est_acf(f_1, 'unbiased')
    acf_unbiased_1 = acf_unbiased_1 / np.max(np.abs(acf_unbiased_1))
   
    ###
    # second function
    ###
    f_2 = np.random.normal(0.0, 1.0, N_2)

    acf_biased_2 = own_spec.est_acf(f_2, 'biased')
    acf_biased_2 = acf_biased_2 / np.max(np.abs(acf_biased_2))  
    
    acf_unbiased_2 = own_spec.est_acf(f_2, 'unbiased')   
    #acf_unbiased_2 = acf_unbiased_2 / np.max(np.abs(acf_unbiased_2)) 
   
    #####
    # estimate psd with periodogram and 2 correlograms from above
    #####   

    # find periodogram by simple fft and abs()**2
    Ome_1 = np.linspace(-np.pi, np.pi, N_freq)
    
    psd_p_1 = own_spec.find_periodogram(f_1, Ome_1)
    psd_p_1 = psd_p_1 / np.max(np.abs(psd_p_1))
    
    psd_c_biased_1 = own_spec.find_correlogram( acf_biased_1, Ome_1 )
    psd_c_biased_1 = psd_c_biased_1 / np.max(np.abs(psd_c_biased_1))
    
    psd_c_unbiased_1 = own_spec.find_correlogram( acf_unbiased_1, Ome_1)
    psd_c_unbiased_1 = psd_c_unbiased_1 / np.max(np.abs(psd_c_unbiased_1))   
   
    # larger N
    Ome_2 = np.linspace(-np.pi, np.pi, len(f_2))

    psd_p_2 = own_spec.find_periodogram(f_2, Ome_2)    
    #psd_p_2 = psd_p_2 / np.max(np.abs(psd_p_2))    
    
    psd_c_biased_2 = own_spec.find_correlogram( acf_biased_2, Ome_2 )
    #psd_c_biased_2 = psd_c_biased_2 / np.max(np.abs(psd_c_biased_2))
    
    psd_c_unbiased_2 = own_spec.find_correlogram( acf_unbiased_2, Ome_2)
    #psd_c_unbiased_2 = psd_c_unbiased_2 / np.max(np.abs(psd_c_unbiased_2))   
    
    # problem: result of psd_unbiased is not a psd (possessing negative values)
    # show why this is true
#    tri = 1-abs(N_acf_2)/N_2
#    tri_reverse = tri[::-1]
#    tri = np.append(tri_reverse[0:len(tri)-1], tri) 
#     
    W = np.zeros((N_1, N_1), dtype=complex) 
    for n in np.arange(0,N_1):
        for m in np.arange(0,N_1):
            W[n,m]=acf_unbiased_1[N_1-1+m-n]

    L, X = np.linalg.eig(W)
    #print np.min(L) 
     
     
    ######################################      
    # plotting
    ######################################    
    font = {'size'   : 26}
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)      
    
    # plot acfs
    plt.figure(1)    
    plt.subplot(121)
    plt.plot(N_acf_1, acf_unbiased_1,'b', label='$L=N-k$')      
    plt.plot(N_acf_1, acf_biased_1,'r', label='$L=N$')
    #plt.plot(N_acf_1, (1-abs(N_acf_1)/N_1)*acf_wgn_unbiased_1,'--')
    plt.xlabel('$k$')
    plt.ylabel('$\hat{r}[k]$')   
    plt.grid(True)    
    plt.title('$N=$'+str(N_1))
    plt.legend(loc='upper right') 
    plt.axis([-N_1, N_1, -.5, 1.1])
    
    plt.subplot(122)
    plt.plot(N_acf_2, acf_unbiased_2,'b', label='$L=N-k$')  
    plt.plot(N_acf_2, acf_biased_2,'r', label='$L=N$')    
    plt.xlabel('$k$')
    plt.ylabel('$\hat{r}[k]$')  
    plt.grid(True)    
    plt.title('$N=$'+str(N_2))    
    plt.legend(loc='upper right') 
    plt.axis([-N_2, N_2, -.5, 1.1])    
    
    # plot psd1
    plt.figure(2)    
    
    plt.subplot(321)    
    plt.plot(Ome_1, psd_p_1)
    #plt.plot(Ome_corr_1, psd_biased_1, label='corr., b')          
    plt.title('$N=$'+str(N_1))    
    plt.grid(True); #plt.legend(loc='upper right') 
    plt.ylabel('$\hat{\Phi}_p(\Omega)$')   
    
    plt.subplot(323)    
    plt.plot(Ome_1, psd_c_biased_1, label='$L=N$')
    plt.grid(True); plt.legend(loc='upper right') 
    plt.ylabel('$\hat{\Phi}_c(\Omega)$')       
    
    plt.subplot(325)        
    plt.plot(Ome_1, psd_c_unbiased_1, label='$L=N-k$')  
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.grid(True); plt.legend(loc='upper right') 
    
    plt.subplot(322)    
    plt.plot(Ome_2, psd_p_2)      
    plt.grid(True); #plt.legend(loc='upper right') 
    plt.title('$N=$'+str(N_2))    
    plt.ylabel('$\hat{\Phi}_p(\Omega)$')   

    plt.subplot(324)    
    plt.plot(Ome_2, psd_c_biased_2, label='$L=N$')      
    plt.grid(True); plt.legend(loc='upper right') 
    plt.ylabel('$\hat{\Phi}_c(\Omega)$')   

    plt.subplot(326)    
    plt.plot(Ome_2, psd_c_unbiased_2, label='$L=N-k$')          
    plt.xlabel('$\Omega$'); plt.ylabel('$\hat{\Phi}_c(\Omega)$')   
    plt.grid(True);  plt.legend(loc='upper right') 
    
    plt.show()


    
########################
# make it executable
########################
if __name__ == "__main__":
    main()