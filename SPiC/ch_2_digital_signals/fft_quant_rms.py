# -*- coding: utf-8 -*-

########################
######## importing
########################
# 
import numpy as np
import matplotlib.pyplot as plt


# Define axes
log_2_N = np.arange(1, 17)   # 2**1 to 2**16
w = np.arange(1, 17)      # 1 to 16

# construct grid
log_2_N_grid, w_grid = np.meshgrid( log_2_N, w)

# sqr of quantization error
def f( log_2_N, w):
    return 1/( 2**log_2_N /12 * 2**(-2.*w) )

Z = f( log_2_N_grid, w_grid)

# Plot
plt.figure(figsize=(10, 5))

plt.subplot(121)
plt.imshow(
    Z,
    aspect='auto',
    origin='lower',
    extent=[log_2_N.min(), log_2_N.max(), w.min(), w.max()]
)

plt.colorbar()#label='RMS')
plt.xlabel('$\\log_2(N)$')
plt.ylabel('$w$')
plt.title('Signal-to-Quantization-Error-Power of an DFT bin (dB)')


plt.subplot(122)
plt.imshow(
    10*np.log10( Z ),
    aspect='auto',
    origin='lower',
    extent=[log_2_N.min(), log_2_N.max(), w.min(), w.max()]
)

plt.colorbar()#label='RMS (dB)')
plt.xlabel('$\\log_2(N)$')
plt.ylabel('$w$')
plt.title('Signal-to-Quantization-Error-Power of an DFT bin (dB)')
plt.tight_layout()
plt.show()