# visualization/spectrum_plot.py

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(sol, N, t_eval):
    """
    Trace le spectre de Fourier des θ_i(t).
    """
    dt = t_eval[1] - t_eval[0]
    freqs = np.fft.rfftfreq(len(t_eval), d=dt)
    
    plt.figure(figsize=(8, 5))
    
    for i in range(N):
        theta_i = sol.y[i]
        fft_theta = np.fft.rfft(theta_i)
        magnitude = np.abs(fft_theta)
        
        plt.plot(freqs, magnitude, label=f"$\\theta_{i+1}$")
    
    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Amplitude")
    plt.title("Spectre des θ_i(t)")
    plt.legend()
    plt.grid()
    plt.show()
