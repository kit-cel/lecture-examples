import numpy as np
import scipy.signal as sig
import scipy.linalg as lin
import matplotlib.pyplot as plt

bpsk_map = np.array([-1, 1])

def generate_random_BPSK_symbols(nsym):
    random_indices = np.random.randint(low=0, high=2, size=nsym)
    return bpsk_map[random_indices]

def detector(rx_sym):
    return np.array([1 if sym > 0 else (-1) for sym in rx_sym])

nsym = 1000000  # number of symbols
signal = generate_random_BPSK_symbols(nsym)

h_c = np.array([1 ,0.7, 0.5, -0.1])

SNR_dB = np.linspace(0,30,50)
sigma2 = 10**(-SNR_dB/10)

BER_n = np.zeros_like(SNR_dB)
BER_zf = np.zeros_like(SNR_dB)
BER_mmse = np.zeros_like(SNR_dB)
i=0
for s in sigma2:
    z = np.sqrt(s) * np.random.randn(signal.size + h_c.size - 1)
    r = np.convolve(h_c, signal) + z
    # no equalization
    r_det = detector(r[:-h_c.size+1])
    BER_n[i] = sum(np.not_equal(r_det, signal)) / r_det.size

    # ZF equalization
    len_e = 3 * h_c.size-1  # number of equalizer coefficients
    len_h = h_c.size  # length of the channel impulse response
    H = lin.toeplitz(np.concatenate((h_c, np.zeros(len_e))), np.concatenate(([h_c[0]], np.zeros(len_e))))
    n = np.zeros(len_h + len_e)
    dirac_pos = 0#len_e//2;
    n[dirac_pos] = 1
    e_ZF = (lin.inv(H.conj().T @ H) @ H.conj().T) @ n
    r_ZF =  np.convolve(e_ZF, r)
    r_det_ZF = detector(r_ZF[dirac_pos: dirac_pos + signal.size])
    BER_zf[i] = sum(np.not_equal(r_det_ZF, signal)) / r_det_ZF.size

    # MMSE equalization
    e_MMSE = (lin.inv(H.conj().T @ H + s * np.eye(H.shape[1])) @ H.conj().T) @ n
    r_MMSE =  np.convolve(e_MMSE, r)
    r_det_MMSE = detector(r_MMSE[dirac_pos: dirac_pos + signal.size])
    BER_mmse[i] = sum(np.not_equal(r_det_MMSE, signal)) / r_det_MMSE.size
    i+=1

plt.figure(dpi=600,figsize=[7,5])
plt.plot(SNR_dB,BER_n,label="No equalization")
plt.plot(SNR_dB,BER_zf, label="Zero forcing")
plt.plot(SNR_dB,BER_mmse, label="MMSE")
plt.legend()
plt.xlabel("SNR (dB)")
plt.ylabel("BER")
plt.yscale('log')
plt.grid()
plt.ylim(1e-4,1)
plt.xlim(0,20)
plt.tight_layout()
plt.savefig("BER_equalization.png")





