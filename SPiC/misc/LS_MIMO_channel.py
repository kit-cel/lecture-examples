# -*- coding: utf-8 -*-
"""
Created on Wed May  7 09:50:25 2014

@author: jaekel
"""


########################
######## importing
########################
  
import numpy as np
import matplotlib.pyplot as plt

import platform
import sys

if platform.system() == "Windows":
    sys.path.append('Z:\\simulations_and_programs\python\comm_sys')
    
if platform.system() == "Linux":    
    sys.path.append('/home/jaekel/INT/programs/python/comm_sys/')
    

import utils_communications as uc


########################
# main function
########################
def main():
    
    # dimension of the mimo system
    N_tx = 2
    N_rx = 2     
    
    # transmitting several valid Tx vectors out of a given constellation and  
    # determining the LS estimation by using the pseudo-inverse of H
    
    # initialize
    M = 4
    mod_scheme = 'QAM'
    
    ser = []    
    ber = []
    constellation = uc.find_constellation(M, mod_scheme)
    
    EbN0_db_min = 0 
    EbN0_db_max = 20
    
    EbN0_db_range = uc.my_range(EbN0_db_min, EbN0_db_max, 1)
    EsN0_db_range = [x+10*np.log10(M) for x in EbN0_db_range]

    max_errors = 1e2
    max_syms = 1e4

   # loop for snr
    for snr_db in EsN0_db_range:
        x_observed = []
        y_observed = []
        num_errors = 0
        num_b_errors = 0
        num_syms = 0
    
        # loop for errors
        while (num_errors<max_errors and num_syms<max_syms):
    
            # data symbols and tx signal
            d = np.random.randint( 0, M, size=N_tx)
            s = uc.modulate(d, constellation) 

            # add fading and noise
            H = 1./np.sqrt(2) * (np.random.normal(size=(N_rx, N_tx)) + 1j*np.random.normal(size=(N_rx, N_tx)) )
            H = H / np.sqrt(np.linalg.norm(H[:,0]))
            H = H / np.linalg.norm(H)
            #H = np.eye(2)

            y_fade = np.dot(H, s)
            noise = uc.add_awgn( y_fade, snr_db, M, 'QAM') - y_fade

            r = y_fade + noise                
                   
            # collecting signals from the first rx antenna            
            y_observed.append(r[0])   

            # demodulate
            x_est = np.dot( np.linalg.pinv(H), r)

            x_observed.append(x_est[0])            
            
            d_est = uc.demod_ML_awgn( x_est, constellation )

            k = int(np.log2(M))

            for n in range(N_tx):
                d_bin = np.binary_repr(d[n], width = k)
                d_est_bin = np.binary_repr(d_est[n], width = k)
                num_b_errors += sum( map( lambda x,y: x!=y, d_bin, d_est_bin) )
                
            # counting errors and symbols
            num_errors += sum( map( lambda x,y: x!=y, d, d_est) )
            num_syms += 1


        # show progress
        done = (snr_db - min(EsN0_db_range))/(max(EsN0_db_range)- min(EsN0_db_range))*100.0
        print 'Done: %3.2f percent' % done            
    
        #print '\t norm(H)=', np.linalg.norm(H), 'P_tx=', np.average(P_tx), 'P_rx=', np.average(P_rx), '\n'    
    
        # ser and ber
        ser.append( num_errors/float(num_syms) )
        ber.append( num_b_errors/(num_syms*np.log2(M)))
        
    
    # theoretical values
    ser_theo = uc.ser_awgn_theory(EbN0_db_range, M, 'QAM')
    

    # plotting
    plt.figure(1)
    plt.plot(EbN0_db_range, ser, 'o', label="simulation SER")
    plt.plot(EbN0_db_range, ber, 'x', label="simulation BER")    
    plt.plot(EbN0_db_range, ser_theo, label="theory AWGN")
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc='upper right') 
    
    plt.xlabel('E_b/N_0 (dB)')
    plt.ylabel('SER')
    
        
    plt.figure(2)
    re=[]
    im=[]
    for c in constellation:
        re.append(c.real)
        im.append(c.imag)

    re_y=[]
    im_y=[]   
    for y in y_observed:
        re_y.append(y.real)
        im_y.append(y.imag)           

    re_x=[]
    im_x=[]   
    for x in x_observed:
        re_x.append(x.real)
        im_x.append(x.imag)                  
           
    plt.plot(im, re, 'ro', im_y, re_y, 'x', im_x, re_x, '+')        
    plt.grid(True)
    plt.xlabel('Inphase')
    plt.ylabel('Quadrature component')        
    plt.show()

########################
# make it executable
########################
if __name__ == "__main__":
    main()
