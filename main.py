import numpy as np
import scipy.signal as sig
from scipy.signal import freqz
import matplotlib.pyplot as plt
# import iir_filter as filter

f_s = 48000    # Sample frequency in Hz
f_c = 2000     # Cut-off frequency in Hz
order = 2    # Order of the butterworth filter

# filter.tim_iirfilter(order, f_c, output='zpk', btype='lowpass', analog=False, fs=f_s, ftype='butter')


sig.iirfilter(order, f_c, output='zpk', btype='lowpass', analog=False, fs=f_s, ftype='butter')

# b, a = sig.butter(order, f_c/24000, analog=False)
#
# sos = sig.tf2sos(b, a, pairing='keep_odd')
#
# w, H = freqz(b, a, 100)                  # Calculate the frequency response
# w *= f_s / (2 * np.pi)                       # Convert from rad/sample to Hz
#
# # Plot the amplitude response
# plt.subplot(2, 1, 1)
# plt.suptitle('Bode Plot')
# H_dB = 20 * np.log10(abs(H))              # Convert modulus of H to dB
# plt.plot(w, H_dB)
# plt.ylabel('Magnitude [dB]')
# plt.xlim(0, f_s / 2)
# plt.ylim(-80, 6)
# plt.axvline(f_c, color='red')
# plt.axhline(-3, linewidth=0.8, color='black', linestyle=':')
# plt.show()
#
# # ===========   Print coefficients   ===========
# print('Coefficients of recursive part:')
# print('A: ')
# print(['%1.10f'%ai for ai in a])
# print('B: ')
# print(['%1.10f'%bi for bi in b])
# print('\n')
# print('Coefficients of sos part:')
# print('Section \t b1 \t\t b2 \t\t b3')
# for n in range(sos.shape[0]):
#     print('%d \t\t %1.10f \t %1.10f \t %1.10f'%(n, sos[n, 0], sos[n, 1], sos[n, 2]))
# print('\n')
# print('Section \t a1 \t\t a2 \t\t a3')
# for n in range(sos.shape[0]):
#     print('%d \t\t %1.10f \t %1.10f \t %1.10f'%(n, sos[n, 3], sos[n, 4], sos[n, 5]))
