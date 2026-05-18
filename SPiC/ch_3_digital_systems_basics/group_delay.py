# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 11:22:06 2015

@author: jaekel
"""


########################
######## importing
########################
from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib

import platform
import sys

if platform.system() == "Windows":
    sys.path.append('Z:\\Simulationen_Programme\kommunikationssysteme')
    
if platform.system() == "Linux":   
    sys.path.append('/home/jaekel/INT/Simulationen_Programme/kommunikationssysteme/')
    sys.path.append('/home/jaekel/Dokumente/Simulationen_Programme/kommunikationssysteme/')
        
import utils_communications as uc 


########################
# main function
########################

# carrier frequency and sampling time
f_c = 2
t_sample = 1/(10*f_c)

#f_c = 0
#t_sample = .1

# parameters for rrc
beta = 0.33
syms_per_filt = 4  # symbols per filter (plus minus in both directions)    

t_symb = 1
n_up = t_symb / t_sample           # samples per symbol such that Nyquist is fulfilled
K_filt = 2*syms_per_filt*n_up+1

# vector of sampled time
t = np.arange(0, K_filt)*t_sample    

# generating rrc pulse 
rrc = uc.get_rrc_ir(K_filt, n_up, t_symb, beta)
#E_rrc = np.sum(np.abs(rrc)**2)
#rrc = rrc/np.sqrt(E_rrc) 

# modulate rrc onto carrier with frequency f_c
carrier = np.cos(2*np.pi*f_c*t)
s_modulated = rrc * carrier

# channels where good guy fulfills necessary symmetrie for linear phase 
# whereas bad guy doesn't 
N_channel = 351
t_extended = np.arange(0, K_filt+N_channel-1)*t_sample    

_h_good = np.random.rand((N_channel-1)//2)
h_good_baseband = np.append( _h_good[::-1], np.append(1, _h_good))
h_good_baseband /= np.sqrt( np.sum( np.abs(h_good_baseband)**2))

s_good_baseband = np.convolve(rrc, h_good_baseband)
s_good= s_good_baseband * np.cos(2 * np.pi*f_c * t_extended)

h_bad_baseband = np.append( np.random.rand((N_channel-1)//2), np.append(1, np.random.rand((N_channel-1)//2) ))
h_bad_baseband /= np.sqrt( np.sum( np.abs(h_bad_baseband)**2))

s_bad_baseband = np.convolve(rrc, h_bad_baseband)
s_bad = s_bad_baseband * np.cos(2 * np.pi*f_c * t_extended)

# frequency regime
H_good = np.fft.fft( np.append(h_good_baseband, np.zeros(5*len(h_good_baseband))))
H_bad = np.fft.fft( np.append(h_bad_baseband, np.zeros(5*len(h_bad_baseband))))

phi_good = np.remainder(np.angle(H_good), 2*np.pi)  - np.pi
phi_bad = np.remainder(np.angle(H_bad), 2*np.pi)  - np.pi

tau_good = np.remainder(np.diff(phi_good), 2*np.pi) - np.pi
tau_bad = np.remainder(np.diff(phi_bad), 2*np.pi) - np.pi

# plotting
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

plt.close('all')

plt.figure(1)
plt.subplot(241)
plt.plot( t, rrc)
plt.plot( t, s_modulated)
plt.grid(True)
plt.title('$s[n]$')

plt.subplot(242)
plt.stem( h_good_baseband)
plt.grid(True)
plt.title('$h[n]$')

plt.subplot(243)
plt.plot( t_extended, s_good)
plt.grid(True)
plt.title('$r[n]$')

plt.subplot(245)
plt.plot( t, rrc)
plt.plot( t, s_modulated)
plt.grid(True)

plt.subplot(246)
plt.stem( h_bad_baseband)
plt.grid(True)

plt.subplot(247)
plt.plot( t_extended, s_bad)
plt.grid(True)
plt.xlabel('$t/s$')

plt.subplot(244)
plt.plot( tau_good )
plt.grid(True)
plt.title('$\\tau_g[k]$')
plt.ylim((-np.pi, np.pi))

plt.subplot(248)
plt.plot( tau_bad ) 
plt.grid(True)
plt.xlabel('$k$')
plt.ylim((-np.pi, np.pi))

plt.show()

# that's it
print('Done!')
