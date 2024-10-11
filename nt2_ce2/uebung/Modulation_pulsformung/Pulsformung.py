# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 10:31:13 2014

NTII Demo - Pulsformung

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

# Ueberabtastung (Samples pro Symbol)
N = 8

# RRC-Filter, Rolloff-Faktor, Anzahl Filterkoeffizienten
alpha = 0.35
N_rrc = N*4+1

# FFT Length
N_FFT = 512

# Pruefe Eingaben
assert (K > 0 and (K & (K - 1)) == 0), 'K muss eine Potenz von 2 sein'
assert (N > 0 and N%2 == 0), 'N muss groesser Null sein und gerade'
assert (alpha >= 0 and alpha <= 1), 'Fuer den Rolloff-Faktor gilt: 0 <= alpha <= 1'
assert (N_rrc > 0 and N_rrc%2 != 0), 'N_rrc muss groesser Null sein und ungerade'


###############################################################################
## Sender
###############################################################################

# QPSK Symbole erzeugen
s = 2 * np.random.randint(2, size=K)-1 + 1j * (2*np.random.randint(2, size=K)-1)

# Leistungsnormalisierung auf 1
s = 1/np.sqrt(2)*s

# Ueberabtasten um Faktor N
s_up = np.zeros(K*N,np.complex)
s_up[::N] = s;

# Rechteckfilter
h_rect = np.ones(N)
s_tx_rect = sig.lfilter(h_rect,1,s_up)
s_tx_rect = 1/np.sqrt(np.mean(np.abs(s_tx_rect)**2))*s_tx_rect
print "Mittlere Leistung von s_tx_rect = %f" % np.mean(np.abs(s_tx_rect)**2)

# Dreiecksfilter
h_tri = np.concatenate([np.arange(0,int(N)+1),np.arange(int(N)-1,0,-1)])
s_tx_tri = sig.lfilter(h_tri,1,s_up)
s_tx_tri = 1/np.sqrt(np.mean(np.abs(s_tx_tri)**2))*s_tx_tri
print "Mittlere Leistung von s_tx_tri  = %f" % np.mean(np.abs(s_tx_tri)**2)

# Root-Raised-Cosine (RRC) Filter 
h_rrc = rrc.get_rrc_ir(N_rrc,N,1.0,alpha)
s_tx_rrc = sig.lfilter(h_rrc,1,s_up)
s_tx_rrc = 1/np.sqrt(np.mean(np.abs(s_tx_rrc)**2))*s_tx_rrc
print "Mittlere Leistung von s_tx_rect = %f" % np.mean(np.abs(s_tx_rrc)**2)


##############################################################################
# Ausgabe
##############################################################################

# Einschwingzeit RRC Filter (Hilfsgroesse)
N_osc = (N_rrc-1)/2


fig1 = plt.figure()
fig1.suptitle("Sender", fontsize=14, fontweight='bold')

ax1 = fig1.add_subplot(2,1,1)
ax1.set_title('Symbole (Tx, Realteil)')
ax1.stem(np.array(range(100)),s.real[:100])
ax1.grid(True)
ax1.set_xlabel('k (t/Ts)')
ax1.set_ylabel('Amplitude')

ax2 = fig1.add_subplot(2,1,2)
ax2.set_title('Symbole ueberabgetastet (Realteil)')
ax2.stem(np.array(range(100)),s_up.real[:100])
ax2.grid(True)
ax2.set_xlabel('k (t/Ts/N)')
ax2.set_ylabel('Amplitude')



fig2 = plt.figure()
fig2.suptitle("Pulsformung (Rechteck)", fontsize=14, fontweight='bold')

ax1 = fig2.add_subplot(1,3,1)
ax1.set_title('Impulsantwort Rechteck')
ax1.stem(np.arange(np.ceil(-N/2),np.ceil(N/2)),h_rect)
ax1.set_xlim(-N_osc,N_osc+1)
ax1.set_ylim(-0.1,1.1)
ax1.grid(True)
ax1.set_xlabel('k (t/Ts/N)')
ax1.set_ylabel('Amplitude')

ax2 = fig2.add_subplot(1,3,2)
ax2.set_title('Signal nach Rechteck-Pulsformung (Realteil, Blau)')
ax2.plot(np.array(range(256)),s_tx_rect.real[0:256])
ax2.plot(np.array(np.arange(0,256,N)),np.real(s[0:np.int(np.ceil(256/N))]),'go',markersize=4)
ax2.set_xlim(0,256)
ax2.grid(True)
ax2.set_xlabel('k (t/Ts/N)')
ax2.set_ylabel('Amplitude')

ax3 = fig2.add_subplot(1,3,3)
ax3.set_title('PSD Rechteck-Pulsformung')
Pxx_rect = 1/(K*N/N_FFT)*np.abs(np.fft.fftshift(np.fft.fft(np.reshape(s_tx_rect,(-1,N_FFT)),axis=1))).sum(0)
f = np.linspace(-0.5,0.5,len(Pxx_rect))
ax3.semilogy(f, Pxx_rect, 'b')
ax3.set_xlim(-0.5,0.5)
ax3.set_ylim([1e-1, np.sqrt(N*K)])
ax3.grid(True)
ax3.set_xlabel('n (f/N/Ts)')

fig3 = plt.figure()
fig3.suptitle("Pulsformung (Dreieck)", fontsize=14, fontweight='bold')

ax1 = fig3.add_subplot(1,3,1)
ax1.set_title('Impulsantwort Dreieck')
ax1.stem(np.arange(np.ceil(-N),np.ceil(N)),h_tri)
ax1.set_xlim(-N_osc,N_osc+1)
ax1.set_ylim(-0.1,1.1*np.max(h_tri))
ax1.grid(True)
ax1.set_xlabel('k (t/Ts/N)')
ax1.set_ylabel('Amplitude')

ax2 = fig3.add_subplot(1,3,2)
ax2.set_title('Signal nach Dreieck-Pulsformung (Realteil, Blau)')
ax2.plot(np.array(range(256)),s_tx_tri.real[N/2:256+N/2])
ax2.plot(np.array(np.arange(0,256,N)),np.real(s[0:np.int(np.ceil(256/N))]),'go',markersize=4)
ax2.set_xlim(0,256)
ax2.grid(True)
ax2.set_xlabel('k (t/Ts/N)')
ax2.set_ylabel('Amplitude')

ax3 = fig3.add_subplot(1,3,3)
ax3.set_title('PSD Dreieck-Pulsformung')
Pxx_tri = 1/(K*N/N_FFT)*np.abs(np.fft.fftshift(np.fft.fft(np.reshape(s_tx_tri,(-1,N_FFT)),axis=1))).sum(0)
f = np.linspace(-0.5,0.5,len(Pxx_tri))
ax3.semilogy(f, Pxx_tri, 'b')
ax3.set_xlim(-0.5,0.5)
ax3.set_ylim([1e-1, np.sqrt(N*K)])
ax3.grid(True)
ax3.set_xlabel('n (f/N/Ts)')



fig4 = plt.figure()
fig4.suptitle("Pulsformung (RRC)", fontsize=14, fontweight='bold')

ax1 = fig4.add_subplot(1,3,1)
ax1.set_title('Impulsantwort RRC')
ax1.stem(np.array(np.arange(-N_osc,N_osc+1)),h_rrc)
ax1.set_xlim(-N_osc,N_osc+1)
ax1.grid(True)
ax1.set_xlabel('k (t/Ts/N)')
ax1.set_ylabel('Amplitude')

ax2 = fig4.add_subplot(1,3,2)
ax2.set_title('Signal nach RRC-Pulsformung (Realteil, Blau)')
ax2.plot(np.array(range(256)),s_tx_rrc.real[(N_rrc-1)/2:256+(N_rrc-1)/2])
ax2.plot(np.array(np.arange(0,256,N)),np.real(s[0:np.int(np.ceil(256/N))]),'go',markersize=4)
ax2.set_xlim(0,256)
ax2.grid(True)
ax2.set_xlabel('k (t/Ts/N)')
ax2.set_ylabel('Amplitude')

ax3 = fig4.add_subplot(1,3,3)
ax3.set_title('PSD RRC-Pulsformung')
Pxx_rrc = 1/(K*N/N_FFT)*np.abs(np.fft.fftshift(np.fft.fft(np.reshape(s_tx_rrc,(-1,N_FFT)),axis=1))).sum(0)
f = np.linspace(-0.5,0.5,len(Pxx_rrc))
ax3.semilogy(f, Pxx_rrc, 'b')
ax3.set_xlim(-0.5,0.5)
ax3.set_ylim([1e-1, np.sqrt(N*K)])
ax3.grid(True)
ax3.set_xlabel('n (f/N/Ts)')


print "Energie von Pxx_rect = %f" % (np.sum(np.abs(Pxx_rect*(f[2]-f[1]))**2))
print "Energie von Pxx_tri  = %f" % (np.sum(np.abs(Pxx_tri*(f[2]-f[1]))**2))
print "Energie von Pxx_rrc  = %f\n" % (np.sum(np.abs(Pxx_rrc*(f[2]-f[1]))**2))

n_start = np.argmin(np.abs(f+0.05))
n_stop  = np.argmin(np.abs(f-0.05))

print "Energie von Pxx_rect (-0.05...+0.05) = %f" % (np.sum(np.abs(Pxx_rect[n_start:n_stop]*(f[2]-f[1]))**2))
print "Energie von Pxx_tri (-0.05...+0.05)  = %f" % (np.sum(np.abs(Pxx_tri[n_start:n_stop]*(f[2]-f[1]))**2))
print "Energie von Pxx_rrc (-0.05...+0.05)  = %f\n" % (np.sum(np.abs(Pxx_rrc[n_start:n_stop]*(f[2]-f[1]))**2))

fig5 = plt.figure()

ax1 = fig5.add_subplot(1,1,1)
ax1.semilogy(f,Pxx_rect, label='Rechteck')
ax1.semilogy(f,Pxx_tri, label='Dreieck')
ax1.semilogy(f,Pxx_rrc, label='RRC')

ax1.grid(True)
ax1.set_xlim(-0.5,0.5)
ax1.set_ylim([1e-1, np.sqrt(N*K)])
ax1.set_xlabel('n (f/N/Ts)')
ax1.set_ylabel('Amplitude')
ax1.set_title('PSD Rechteck-/Dreieck/RRC-Pulsformung')
ax1.legend()

plt.show()
