import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

c_output = pd.read_csv("/home/ogutierreze/GPDS_Proyects/ANE2/ANE2-GCPDS/Sw/Software-HackRF/ANE2-driver/Samples/output2.csv")

f_c = np.array(c_output['Frequency_Hz'].values)
Pxx_c = np.array(c_output['PSD'].values)

plt.figure(figsize=(20, 10))
plt.semilogy(f_c, Pxx_c)
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (dB/Hz)')
plt.title('Power Spectral Density (C Output)')
plt.axhline(y=10**(-7.5), color='r', linestyle='--', linewidth=2, label='Línea horizontal')
plt.grid()
plt.show()
