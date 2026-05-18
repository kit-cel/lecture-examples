# -*- coding: utf-8 -*-
"""
Created on Thu May 07 09:22:52 2015

@author: jaekel
"""


########################
######## importing
########################
  
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import platform
import sys

if platform.system() == "Windows":
    sys.path.append('Z:\\Simulationen_Programme\Kommunikationssysteme')
    
if platform.system() == "Linux":   
    sys.path.append('/home/jaekel/INT/Simulationen_Programme/Kommunikationssysteme/')
    sys.path.append('/home/jaekel/Dokumente/Simulationen_Programme/Kommunikationssysteme/')
        
import utils_communications as uc  

plt.close('all')

########################
# main function
########################
   
# define channel impulse response
while True:
    N_h = 15 #np.random.randint(25)
    if N_h>0:
        break

h = np.random.randn(N_h) + 1j * np.random.randn(N_h)
h /= np.linalg.norm(h)

    
# modulation scheme and constellation points
M = 2
mod_scheme = 'ASK'
constellation = uc.find_constellation(M, mod_scheme)


# generate PN sequence
pn_taps = [1, 0, 1, 0, 0, 1]      # mind that MSB (corresponding to the highest exponent) comes last
pn_seed = [0, 0, 0, 0, 1] 

pn = uc.lfsr(pn_taps, pn_seed)  
pn_mod = uc.modulate( pn, constellation ) 


# define transmission matrix X in the model y=Xh+n
X = np.zeros( (len(pn)+N_h-1, N_h) )

for k in np.arange(N_h):        
    X[k:k+len(pn), k] = pn_mod


# filter signal and add noise 
snr_db = np.array([0, 10, 20])

# initialze result of estimation
h_est = np.zeros( (len(snr_db), N_h), dtype=complex )

for k, s in enumerate(snr_db):
    x = np.convolve(pn_mod, h)
    y = uc.add_awgn( x, s, M, mod_scheme)

    # do channel estimation by LS estimation
    h_est[ k, : ] = np.dot( np.linalg.pinv(X), y )



# plotting
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)  

plt.figure(2)    
for k, s in enumerate(snr_db):
    plt.subplot(1, len(snr_db), k+1 )
    plt.plot( np.real(h), label='$h[n]$' ) 
    plt.plot( np.real(h_est[ k, :]), label='$h_{\mathrm{est}}[n]$') 
    plt.grid(True)
    plt.legend( loc='upper right' )
    plt.xlabel('$n$')
    plt.title('SNR = '+str(s)+' (dB)')


plt.show()


# that's it folks!
print('Done!')

