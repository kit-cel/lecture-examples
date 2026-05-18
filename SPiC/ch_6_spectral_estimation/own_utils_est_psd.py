# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 09:13:45 2014

@author: jaekel
"""

import numpy as np



########################
# own acf estimator
########################
def est_acf(y, est_type):
    """
    estimates acf given a number of observation
    
    Remark: signal is assumed to be starting from 0 to length(y)-1
    
    IN: observations y, est_type (biased / unbiased)
    OUT: estimated acf, centered around 0
    """
    N = np.size(y)
    r = np.zeros(N, dtype=complex)
    
    # loop lags of acf
    for k in np.arange(0,N):
        temp = np.sum( y[k:N] * np.conjugate(y[0:(N-k)]) )

        # type of estimator
        if est_type == 'biased':
            r[k] = temp/N
        elif est_type == 'unbiased':
            r[k] = temp/(N-k)
        
    # find values for negative indices
    r_reverse = np.conjugate(r[::-1])     
   
    return  np.append(r_reverse[0:len(r)-1], r)  


def est_acf(y, est_type):
    """
    estimates acf given a number of observation
    
    Remark: signal is assumed to be starting from 0 to length(y)-1
    
    IN: observations y, est_type (biased / unbiased)
    OUT: estimated acf, centered around 0
    """
    N = np.size(y)
    r = np.zeros(N, dtype=complex)
    
    # loop lags of acf
    for k in np.arange(0,N):
        temp = np.sum( y[k:N] * np.conjugate(y[0:(N-k)]) )

        # type of estimator
        if est_type == 'biased':
            r[k] = temp/N
        elif est_type == 'unbiased':
            r[k] = temp/(N-k)
        
    # find values for negative indices
    r_reverse = np.conjugate(r[::-1])     
   
    return  np.append(r_reverse[0:len(r)-1], r)  


########################
# own correlogram estimator
########################
def find_correlogram(r, omega):
    """
    estimates correlogram out of the given acf at the frequencies specified in omega
    
    Remark: acf ios assumed to be centered around 0
    
    IN: acf r, frequencies
    OUT: psd
    """
    corr = np.zeros(len(omega), dtype=complex)
  
    N=(len(r)+1)/2
        
    for p in np.arange(-(N-1), (N-1)+1):
        corr += r[p+(N-1)]*np.exp(-1j*omega*p)

    #return np.abs(corr)
    return corr


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
# Bartlett periodogram estimator
########################
def find_bartlett_estimate(y, M, omega):
    """
    estimates periodogram out of the given observation at the frequencies specified in omega
    using Bartlett's method
    
    IN: observation y, group size M, frequencies Omega
    OUT: psd estimator
    """
    
    N = len(y)
    K = int( float(N)/M )
    
    per = np.zeros(len(omega), dtype=complex)
    
    k = 0
    while k<K:
        
        yk = y[ k*M : (k+1)*M ]         # mind that the upper limit is not included    
        Yk = find_periodogram(yk, omega)
        per = 1.0/(k+1) * ( k*per + Yk )  

        k += 1            

    return per
    
    
########################
# Welch periodogram estimator
########################
def find_welch_estimate(y, M, O, omega, window=[]):
    """
    estimates periodogram out of the given observation at the frequencies specified in omega
    using Welch's method
    
    IN: observation y, group size M, overlap O, frequencies Omega, 
        window of length M chosen as rect if not specified otherwise
    OUT: psd estimator
    """
    
    N = len(y)
    K = int( N / (M-O) )
    
    per = np.zeros(len(omega))
    
    if window==[]:
        window = np.ones(M)
    
    window = window / (1.0/M * np.sum(window**2))

    # loop for segments    
    k = 0
    while k<K-1:        
        yk = y[ k*(M-O) : (k*(M-O)+M) ] # mind that the upper limit is not included        
        yk = yk * window        
    

        
        Yk = find_periodogram(yk, omega)

        per = 1.0/(k+1) * ( k*per + Yk )  

        k += 1            

    return per    
    
    

########################
# own Yule Walker periodogram estimator
########################
def find_yulewalker(y, order, s2, omega):
    """
    estimates periodogram using the Yule-Walker method
    
    IN: acf r, YW oder q, noise variance, frequencies
    OUT: psd estimator
    """

    per = np.zeros(len(omega), dtype=complex)

    r = est_acf(y, 'biased') 
    #r = r/np.max(acf)
    N = int( (len(r)+1)/2 )

    # get matrix R for Yule-Walker     
    # note that R is not the autocorrelation matrix, but R = (ACF matrix)^*
    R = np.zeros([order+1, order+1], dtype=complex)        
    for p in range(0, order+1):
        R[:,p] = r[N-1-p : N-1-p+order+1 ]


    # find and solve linear equation system for the coefficients
    b = np.matrix(np.append(s2, np.zeros(order))).T
    theta = np.linalg.solve(R, b)
    theta = np.array(theta / theta[0])

    # get frequency response    
    for p in np.arange(0, order+1):
        per += theta[p]*np.exp(-1j*omega*p)
        
    return s2 / abs(per)**2


    
    
########################
# MUSIC periodogram estimator
########################
def find_music_estimate(y, n, m, omega):
    """
    estimates periodogram out of the given acf at the frequencies specified in omega
    using MUSIC method
    
    IN: signal y, order n, observation number m, frequencies Omega
    OUT: psd estimator
    """
    
    N = len(y)
    per = np.zeros(len(omega), dtype=complex)

    # estimating the autocorrelation matrix

                # THIS realization is a little bit slower
                #    R1 = np.zeros([m,m], dtype=complex)
                #    for t in np.arange(m, N):
                #        yt = (y[(t-m+1):t+1])[::-1]
                #        R1 += np.outer(yt, yt.conjugate() )
                #    R1 = R1/N

    
    # estimating the autocorrelation matrix    
    R = np.zeros([m,m], dtype=complex)
    r = est_acf(y, 'biased')
    for p in np.arange(0, m):
        R[p,:] = r[N-1-p : N-1-p+m ]
     
    per = np.zeros(len(omega), dtype=complex)

    # find eigendecomposition and assign noise space matrix G    
    Lambda, X = np.linalg.eig(R)
    
    indices = Lambda.argsort()[0:m-n]
    G= X[:, indices]
    
    # find pseudospectrum
        # construcing matrix A, containing the vectors a(omega) for
        # samples in the frequency regime given by omega
    A = np.ones([m, len(omega)], dtype=complex)        
    row_freq = np.exp(-1j*omega)
    for p in np.arange(1, m):
        A[p, :] = A[p-1, :] * row_freq

        # matrix product 
    GA = np.array( np.matrix(G).H * np.matrix(A))
    ps = np.sum( np.abs(GA)**2, axis=0 )
    
    per = 1.0 / ps
     
    return per
        
        
        
   
########################
# ESPRIT periodogram estimator
########################
def find_esprit_estimate(y, n, m, omega):
    """
    estimates periodogram out of the given acf at the frequencies specified in omega
    using ESPRIT method
    
    IN: signal y, order n, observation number m, frequencies Omega
    OUT: psd estimator
    """
    
    N = len(y)
    per = np.ones(len(omega), dtype=complex)

    # estimating the autocorrelation matrix    
    R = np.zeros([m,m], dtype=complex)
    r = est_acf(y, 'biased')
    for p in np.arange(0, m):
        R[p,:] = r[N-1-p : N-1-p+m ]
   
    # find eigendecomposition and assign signal space matrix S
            # second operation reverses order to have largest EV first     
    Lambda_R, X_R = np.linalg.eig(R)
    indices = Lambda_R.argsort()[m-n:m][::-1]    
    S= X_R[:, indices]

    # eliminating first and last row
    shape = np.shape(S)
    S1 = np.matrix( S[ :shape[0]-1, :] )
    S2 = np.matrix( S[ 1: , :] )
    
    # constructing matrix Phi = (S1^* S1)^-1 S1^* S2
    # and finding eigenvalues
    Phi = np.linalg.pinv( S1 ) * S2  
    Lambda_Phi, X_Phi = np.linalg.eig(Phi)
  
    # frequencies and amplitudes as angle and abs of the eigenvalues  
    freq_est = -np.angle(Lambda_Phi)
    ampl_est = np.abs(Lambda_Phi)
      
    # find denominator by multiplying contribution of each eigenvalue
    for p in np.arange(0,n):
        per *= ( np.exp(1j*omega) - ampl_est[p]*np.exp(1j*freq_est[p]) )

    return 1/(np.abs(per)**2)
        