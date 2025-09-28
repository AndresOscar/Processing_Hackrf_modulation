import math
import cmath
import struct
import matplotlib.pyplot as plt
import numpy as np

"""--------------------Functions--------------------"""

def fft_scratch(x):
    """
    A custom, from-scratch implementation of the Cooley-Tukey FFT algorithm.
    It works for input arrays whose length is a power of 2.
    """
    N = len(x)
    if N <= 1:
        return x

    # Check if N is a power of 2
    if (N & (N - 1)) != 0:
        raise ValueError("FFT from scratch only works for signal lengths that are powers of 2.")

    even = fft_scratch(x[0::2])
    odd = fft_scratch(x[1::2])

    T = [0] * (N // 2)
    for k in range(N // 2):
        # Calculate twiddle factor
        twiddle = cmath.exp(-2j * cmath.pi * k / N)
        T[k] = twiddle * odd[k]

    X = [0] * N
    for k in range(N // 2):
        X[k] = even[k] + T[k]
        X[k + N // 2] = even[k] - T[k]

    return X

def generate_window(window_type, size, beta=None, alpha=None):
    """
    Generates a specified window function of a given size.
    """
    if size <= 0:
        raise ValueError("Window size must be a positive integer.")

    window = [0.0] * size

    if window_type == 'hamming':
        for n in range(size):
            window[n] = 0.54 - 0.46 * math.cos(2 * math.pi * n / (size - 1))
    elif window_type == 'hann':
        for n in range(size):
            window[n] = 0.5 * (1 - math.cos(2 * math.pi * n / (size - 1)))
    elif window_type == 'rectangular':
        window = [1.0] * size
    elif window_type == 'blackman':
        for n in range(size):
            window[n] = 0.42 - 0.5 * math.cos(2 * math.pi * n / (size - 1)) + 0.08 * math.cos(4 * math.pi * n / (size - 1))
    elif window_type == 'flat_top':
        a0, a1, a2, a3, a4 = 0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368
        for n in range(size):
            window[n] = a0 - a1 * math.cos(2 * math.pi * n / (size - 1)) + a2 * math.cos(4 * math.pi * n / (size - 1)) - a3 * math.cos(6 * math.pi * n / (size - 1)) + a4 * math.cos(8 * math.pi * n / (size - 1))
    elif window_type == 'kaiser':
        if beta is None:
            raise ValueError("Kaiser window requires a 'beta' parameter.")
        I_0_beta = I_0(beta)
        for n in range(size):
            arg = beta * math.sqrt(1 - (2.0 * n / (size - 1) - 1.0)**2)
            window[n] = I_0(arg) / I_0_beta
    elif window_type == 'tukey':
        if alpha is None or not 0 <= alpha <= 1:
            raise ValueError("Tukey window requires an 'alpha' parameter between 0 and 1.")
        half_alpha_N = alpha * (size - 1) / 2.0
        for n in range(size):
            if n < half_alpha_N:
                window[n] = 0.5 * (1 + math.cos(math.pi * (2 * n / (alpha * (size - 1)) - 1)))
            elif n > (size - 1) - half_alpha_N:
                window[n] = 0.5 * (1 + math.cos(math.pi * (2 * n / (alpha * (size - 1)) - 2 / alpha + 1)))
            else:
                window[n] = 1.0
    elif window_type == 'bartlett':
        for n in range(size):
            window[n] = 1 - abs(2 * n / (size - 1) - 1)
    else:
        raise ValueError(f"Unsupported window type: {window_type}")

    return window

def I_0(x):
    """
    Computes the modified zeroth-order Bessel function of the first kind.
    """
    sum_val = 1.0
    term = 1.0
    k = 1
    while term > 1e-9 * sum_val:
        term = term * (x * x / 4.0) / (k * k)
        sum_val += term
        k += 1
    return sum_val

def welch_psd_scratch(signal, fs=1.0, segment_length=256, overlap=0.5, window_type='hamming', **kwargs):
    """
    Calculates the Welch Power Spectral Density (PSD) from scratch.
    """
    if (segment_length & (segment_length - 1)) != 0:
        raise ValueError("For this scratch implementation, segment_length must be a power of 2.")

    N = len(signal)
    step = int(segment_length * (1 - overlap))

    window = generate_window(window_type, segment_length, **kwargs)
    U = sum(w ** 2 for w in window)
    K = (N - segment_length) // step + 1
    P_welch = [0.0] * segment_length

    for k in range(K):
        start = k * step
        segment = signal[start:start + segment_length]
        windowed_segment = [segment[i] * window[i] for i in range(segment_length)]
        X_k = fft_scratch(windowed_segment)
        P_k = [(1 / (fs * U)) * abs(x)**2 for x in X_k]

        for i in range(segment_length):
            P_welch[i] += P_k[i]

    P_welch = [p / K for p in P_welch]

    f = [0.0] * segment_length
    for i in range(segment_length):
        if i < segment_length / 2:
            f[i] = i * fs / segment_length
        else:
            f[i] = (i - segment_length) * fs / segment_length

    return f, P_welch

def moving_average(x, window):
    """Moving average filter"""
    if window <= 1:
        return x[:]
    N = len(x)
    out = [0.0]*N
    half = window//2
    for i in range(N):
        l = max(0, i-half)
        r = min(N-1, i+half)
        s = 0.0
        count = 0
        for j in range(l, r+1):
            s += x[j]
            count += 1
        out[i] = s / max(1, count)
    return out

def moving_median(x, window):
    """Moving median filter"""
    if window <= 1:
        return x[:]
    N = len(x)
    out = [0.0]*N
    half = window//2
    for i in range(N):
        l = max(0, i-half)
        r = min(N-1, i+half)
        seg = x[l:r+1]
        seg_sorted = sorted(seg)
        m = seg_sorted[len(seg_sorted)//2]
        out[i] = m
    return out

def parabolic_interpolate_bin(y_m1, y0, y_p1):
    """Parabolic interpolation for peak refinement"""
    denom = (y_m1 - 2*y0 + y_p1)
    if denom == 0:
        return 0.0
    return 0.5 * (y_m1 - y_p1) / denom

def find_peaks_improved(f, P, smooth_win=5, noise_win=501, thresh_db=6.0, 
                       min_distance_hz=0.0, min_peak_width_bins=1):
    """Improved peak detection algorithm"""
    N = len(P)
    if N < 3:
        return []

    f = list(f)
    P = list(P)
    df = abs(f[1] - f[0]) if N > 1 else 1.0
    min_dist_bins = int(max(1, round(min_distance_hz / df))) if min_distance_hz>0 else 1

    # Smoothing and noise floor calculation
    P_smooth = moving_average(P, smooth_win)
    noise_floor = moving_median(P_smooth, noise_win)
    threshold_linear = [nf * (10**(thresh_db / 10.0)) for nf in noise_floor]

    # Peak detection
    candidate_idxs = []
    for i in range(1, N-1):
        if P_smooth[i] > P_smooth[i-1] and P_smooth[i] > P_smooth[i+1]:
            if P_smooth[i] > threshold_linear[i]:
                candidate_idxs.append(i)

    # Filter nearby peaks
    if candidate_idxs:
        candidate_idxs.sort()
        filtered = []
        for idx in candidate_idxs:
            if filtered and (idx - filtered[-1]) < min_dist_bins:
                prev = filtered[-1]
                if P_smooth[idx] > P_smooth[prev]:
                    filtered[-1] = idx
            else:
                filtered.append(idx)
        candidate_idxs = filtered

    # Refine peaks
    peaks = []
    for i in candidate_idxs:
        if i-1 < 0 or i+1 >= N:
            continue
        delta = parabolic_interpolate_bin(P_smooth[i-1], P_smooth[i], P_smooth[i+1])
        f_peak = f[i] + delta * df
        local_floor = noise_floor[i]
        P_peak = P_smooth[i]
        snr_db = 10*math.log10(P_peak / local_floor) if local_floor>0 else float('inf')
        peaks.append({'idx': i, 'freq': f_peak, 'power': P_peak, 'snr_db': snr_db})

    peaks.sort(key=lambda x: x['freq'])
    return peaks, noise_floor, P_smooth

def read_cs8_file(file_path, num_samples=None):
    """
    Reads a CS8 file (8-bit complex samples: I,Q,I,Q,...)
    Returns a list of complex numbers
    """
    signal = []
    with open(file_path, 'rb') as f:
        data = f.read()
        if num_samples is not None:
            data = data[:num_samples * 2]  # 2 bytes per complex sample
        
        # Convert to complex numbers
        for i in range(0, len(data), 2):
            if i + 1 < len(data):
                # Convert from signed 8-bit to float (-128 to 127 -> -1.0 to 1.0)
                i_val = (data[i] - 128) / 128.0 if data[i] > 127 else data[i] / 128.0
                q_val = (data[i+1] - 128) / 128.0 if data[i+1] > 127 else data[i+1] / 128.0
                signal.append(complex(i_val, q_val))
    
    return signal

# Constants
FS = 20e6
NPERSEG = 4096
OVERLAP = 0.5
WINDOW = 'hamming'

# Ruta específica al archivo CS8
CS8_FILE_PATH = "/home/ogutierreze/GPDS_Proyects/ANE2/ANE2-GCPDS/Sw/Software-HackRF/ANE2-driver/Samples/0.cs8"  # Cambia esta ruta

"""-----------------------------Carga de la señal desde archivo CS8----------------------------------"""

try:
    # Cargar señal desde archivo CS8 (puedes ajustar num_samples si es muy grande)
    signal = read_cs8_file(CS8_FILE_PATH, num_samples=100000)  # Limitar a 100k muestras si es necesario
    print(f"Señal cargada correctamente con {len(signal)} muestras.")
    
except FileNotFoundError:
    print(f"Error: No se encontró el archivo {CS8_FILE_PATH}")
    exit()
except Exception as e:
    print(f"Error al leer el archivo: {e}")
    exit()

"""---------------------------------Procesamiento Principal---------------------------------------"""

# Calcular PSD
f, P_welch = welch_psd_scratch(signal, fs=FS, segment_length=NPERSEG, overlap=OVERLAP, window_type=WINDOW)

# Plot PSD
plt.figure(figsize=(12,6))
plt.semilogy(f, P_welch)
plt.xlabel('Frequency (Hz)', fontsize=12)
plt.ylabel('Power (dB)', fontsize=12)
plt.title("Welch PSD", fontsize=14)
plt.grid(True, which="both", linestyle="--", linewidth=0.7)
plt.show()

# Detectar picos
peaks, noise_floor, P_smooth = find_peaks_improved(
    f, P_welch,
    smooth_win=7,
    noise_win=2001,
    thresh_db=8.0,
    min_distance_hz=1000,
    min_peak_width_bins=2
)

print(f"Detectados {len(peaks)} picos:")
for p in peaks:
    print(f" - {p['freq']:.2f} Hz, power={p['power']:.3e}, SNR={p['snr_db']:.2f} dB")

# Plot con detección de picos
plt.figure(figsize=(14,6))
plt.semilogy(f, P_welch, label='PSD raw')
plt.semilogy(f, P_smooth, alpha=0.9, label='PSD suavizada')

# Umbral de ruido
noise_floor_dB = [10*math.log10(nf) if nf>0 else -300 for nf in noise_floor]
threshold_dB = [val+1 for val in noise_floor_dB]
plt.plot(f, [10**(tdb/10.0) for tdb in threshold_dB], 'r--', linewidth=1.5, label='Umbral ruido')

plt.xlabel('Frecuencia (Hz)')
plt.ylabel('PSD (lineal)')
plt.title('Welch PSD - detección de picos con ancho de banda')
plt.grid(True, which='both', linestyle='--', linewidth=0.6)

# Marcar picos y calcular anchos de banda
for p in peaks:
    idx = p['idx']
    plt.plot(p['freq'], p['power'], marker='o', markersize=7,
             markeredgecolor='k', markerfacecolor='none')
    plt.text(p['freq'], p['power']*1.3, f"{p['freq']:.1f} Hz\n{p['snr_db']:.1f} dB",
             ha='center', va='bottom', fontsize=9)

    # Calcular ancho de banda
    left = idx
    while left > 0 and P_smooth[left] > noise_floor[left]:
        left -= 1
    right = idx
    while right < len(P_smooth)-1 and P_smooth[right] > noise_floor[right]:
        right += 1

    f_left = f[left]
    f_right = f[right]
    bw = f_right - f_left

    print(f"   BW ~ {bw:.2f} Hz (desde {f_left:.2f} Hz hasta {f_right:.2f} Hz)")

    # Graficar ancho de banda
    plt.axvline(f_left, color="gray", linestyle="--", alpha=0.7)
    plt.axvline(f_right, color="gray", linestyle="--", alpha=0.7)
    plt.fill_between(f[left:right], P_smooth[left:right], noise_floor[left:right],
                     color="yellow", alpha=0.3)

plt.legend()
plt.show()