# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 08:32:43 2015

@author: jaekel
"""

########################
######## importing
########################
  
from scipy import signal

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


########################
# main function
########################

N = 32

x = np.ones(N)

conv_x = np.convolve(x, x, 'full')
CONV_x = np.fft.fft(conv_x)

CONV_x_cyclic = np.fft.fft(x)**2
conv_x_cyclic = np.fft.ifft( CONV_x_cyclic )

f_dis = np.arange(-len(conv_x)/2, len(conv_x)/2)
f_dis_cyclic = np.arange(-len(conv_x_cyclic)/2, len(conv_x_cyclic)/2)


# plotting
font = {'size'   : 26}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)

   

plt.figure(1)
plt.subplot(221)
plt.stem( range(np.size(conv_x)), conv_x, label='$x[n]\\ast x[n]$')
plt.grid(True)
plt.xlabel('$n$')
#plt.ylabel('$x[n]*x[n]$')
plt.title('azyklisch')
plt.legend(loc='upper right')

plt.subplot(222)
plt.stem( range(np.size(conv_x_cyclic)), conv_x_cyclic, label='$IFFT\{X[k]\cdot X[k]\}$')
plt.grid(True)
plt.xlabel('$n$')
#plt.ylabel('$IFFT\{ X[k] \cdot X[k] \}$')
plt.title('zyklisch')
plt.legend(loc='upper right')

plt.subplot(223)
plt.stem( f_dis, np.abs( np.fft.fftshift( CONV_x ) ), label='$| FFT\{x[n]*x[n]\}|$' )
plt.grid(True)
plt.xlabel('$k$')
#plt.ylabel('$| FFT\{x[n]*x[n]\}|$')
plt.legend(loc='upper right')


plt.subplot(224)
plt.stem( f_dis_cyclic, np.abs( np.fft.fftshift( CONV_x_cyclic ) ), label='$| X[k]\cdot X[k]|$' )
plt.grid(True)
plt.xlabel('$k$')
#plt.ylabel('$|X[k]\cdot X[k]|$')
plt.legend(loc='upper right')

plt.show()


############
print('Done!')