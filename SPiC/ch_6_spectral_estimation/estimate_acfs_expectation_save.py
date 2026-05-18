# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 12:56:08 2014

@author: jaekel
"""

########################
# illustrating asymptotically unbiased estimation of psd
########################

# importing relevant stuff
import numpy as np
import own_utils_est_psd as own_spec

########################
# main function
########################
def main():
    # parameters: number of points in first resp. second function and in frequency domain
    N = 1e1
    N_acf = np.arange(-N+1, N, 1)
    
    # number of freq. points and freq. range
    N_freq = 16
    Ome = np.linspace(-np.pi, np.pi, N_freq)
    
    # number of realizations for averaging    
    N_aver = 1e2
    
    ###
    # first function
    ###
    # relative width of rectangular 
    relative_width_rect = 0.3
    N_rect = int(N*relative_width_rect)
    f = np.append( np.ones( N_rect), np.zeros( N-N_rect))
    
    
    # loop for realizations
    acf_biased = np.zeros(len(N_acf))
    acf_unbiased = np.zeros(len(N_acf))
    
    psd_p = np.zeros(len(Ome))
    psd_c_biased = np.zeros(len(Ome))
    psd_c_unbiased = np.zeros(len(Ome))
    
    n=0
    while n < N_aver:
    
        # first function    
        #f_1 = f + np.random.normal(0.0, 0.2, N_1)
        f_1 = np.random.normal(0.0, 1.0, N)
        #f_1 = np.sin(1.5*np.arange(0, N_1)) + np.random.normal(0.0, 0.2, N_1)
        
        acf_biased = 1./(n+1) *(float(n)*acf_biased + own_spec.est_acf(f_1, 'biased') )
        acf_unbiased = 1./(n+1) *(float(n)*acf_unbiased + own_spec.est_acf(f_1, 'unbiased') )
        
     
        # find periodogram by simple fft and abs()**2
        psd_p = 1./(n+1) *(n*psd_p + own_spec.find_periodogram(f_1, Ome) )
        
        psd_c_biased = 1./(n+1) *(n*psd_c_biased + own_spec.find_correlogram( acf_biased, Ome) )
        psd_c_unbiased = 1./(n+1) *(n*psd_c_unbiased + own_spec.find_correlogram( acf_unbiased, Ome) )
   
        n += 1
    
        # show progress
        done = float(n)/N_aver*100.0
        print 'Done: %3.2f percent' % done
     
     
    # write to file
    file_name = 'results'

    #file_obj = open( file_name, 'w')
    #file_obj.write('- acf biased - acf unbiased - psd perdiodogram - psd correlogram biased - psd correlogram unbiased \n')
    
#    file_obj.write(str(N_acf) + '; ' + str(acf_biased)+';;;')
#    file_obj.write(str(N_acf) + '; ' +  str(acf_unbiased)+';;;')    
#    file_obj.write(str(Ome) + '; ' +  str(psd_p)+';;;')        
#    file_obj.write(str(Ome) + '; ' +  str(psd_c_biased)+';;;')            
#    file_obj.write(str(Ome) + '; ' +  str(psd_c_unbiased)+';;;')                
#    file_obj.close()     
    
    np.save(file_name, N_acf)
    np.save(file_name, acf_biased)    


    
    

        


    
########################
# make it executable
########################
if __name__ == "__main__":
    main()