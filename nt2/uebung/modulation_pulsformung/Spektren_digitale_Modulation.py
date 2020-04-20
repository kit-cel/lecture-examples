# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 10:31:13 2014

NTII Demo - Quellencodierung - Auswirkungen auf Spektrum des Sendesignals

Systemmodell: Quelle --> QPSK --> Pulsformung

@author: Michael Schwall
"""

from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import scipy.signal as sig
import rrc as rrc
plt.close("all")


###############################################################################
## Systemparameter
###############################################################################

# Anzahl der simulierten Symbole
K = 65536

# Wahrscheinlichkeit fuer ein 1-Bit
P_b_1 = np.array([0.5,0.1])

# Ueberabtastung (Samples pro Symbol)
N = 4

# RRC-Filter, Rolloff-Faktor, Anzahl Filterkoeffizienten
alpha = 0
N_rrc = N*16+1

# FFT Length
N_FFT = 1024

# Pruefe Eingaben
assert (K > 0 and (K & (K - 1)) == 0), 'K muss eine Potenz von 2 sein'
assert (N > 0 and N%2 == 0), 'N muss groesser Null sein und gerade'
assert (alpha >= 0 and alpha <= 1), 'Fuer den Rolloff-Faktor gilt: 0 <= alpha <= 1'
assert (N_rrc > 0 and N_rrc%2 != 0), 'N_rrc muss groesser Null sein und ungerade'


###############################################################################
## Sender
###############################################################################

idx=0
s_tx_rrc = np.zeros((K*N,len(P_b_1)))

while idx < len(P_b_1):
    # Bits erzeugen
    b = (P_b_1[idx]+np.random.uniform(-1.0,0.0,size=K) >= 0).astype(int)
        
    # BPSK Symbole erzeugen
    I = (2*b-1)
        
    print "P(b=1)=%0.2f --> E{I} = %0.2f --> Var{I} = %0.2f" % (P_b_1[idx], I.mean(), I.var())
    
    # Ueberabtasten um Faktor N
    s_up = np.zeros(K*N)
    s_up[::N] = I;
    
    # Root-Raised-Cosine (RRC) Filter 
    h_rrc = rrc.get_rrc_ir(N_rrc,N,1.0,alpha)
    s_tx_rrc[:,idx] = sig.lfilter(h_rrc,1.0,s_up)
    idx += 1
    


##############################################################################
# Ausgabe
##############################################################################

# Einschwingzeit RRC Filter (Hilfsgroesse)
N_osc = (N_rrc-1)/2


fig1 = plt.figure()
fig1.suptitle("Pulsformung (RRC, alpha=%0.2f)" % alpha, fontsize=14, fontweight='bold')

ax1 = fig1.add_subplot(1,2,1)
ax1.set_title('Impulsantwort RRC')
ax1.stem(np.array(np.arange(-N_osc,N_osc+1)),h_rrc)
ax1.set_xlim(-N_osc,N_osc+1)
ax1.grid(True)
ax1.set_xlabel('k (t/Ts/N)')
ax1.set_ylabel('Amplitude')

ax2 = fig1.add_subplot(1,2,2)
ax2.set_title('PSD QPSK mit RRC-Pulsformung')

idx=0
while idx < len(P_b_1):
    Pxx_rrc = 1/N_FFT*(np.abs(np.fft.fftshift(np.fft.fft(np.reshape(s_tx_rrc[:,idx],(-1,N_FFT)),axis=1)))**2).sum(0)
    f = np.linspace(-0.5,0.5,len(Pxx_rrc))
    ax2.plot(f, 10*np.log10(Pxx_rrc))
    idx += 1
    
start, end = ax2.get_ylim()
ax2.yaxis.set_ticks(np.arange(start, end, 10))
    
ax2.set_xlim(-0.5,0.5)
ax2.grid(True)
ax2.set_xlabel('n (f/N/Ts)')
ax2.set_ylabel('Amplitude [dB]')

plt.show()
