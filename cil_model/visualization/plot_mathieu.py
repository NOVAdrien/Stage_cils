# visualization/plot_mathieu.py

import matplotlib.pyplot as plt
import numpy as np
from models.passive_chain_friction_with_activation import PassiveChainFrictionWithActivation

def plot_mathieu(masses, lengths, widths, stiffnesses, k_spatial, phi, alpha, gamma, initial_conditions, t_span, N, t_eval) :

    # Paramètres du colorbar à modifier si besoin pour centrer sur la langue de Mathieu
    As = np.linspace(0.05, 0.3, 4)
    omegas = np.linspace(0.1, 5, 4)
    amplitudes_max = np.zeros((len(As), len(omegas)))

    for i, A_val in enumerate(As):
        A_val = [A_val] * N
        for j, omega_val in enumerate(omegas):
            omega_val = [omega_val] * N
            chain_friction_activation = PassiveChainFrictionWithActivation(
                N, masses, lengths, widths, stiffnesses, A_val, omega_val, k_spatial, phi,
                alpha, gamma,
                initial_conditions
            )

            sol = chain_friction_activation.solve(t_span, t_eval=t_eval, method='RK45')
            amplitudes_max[i, j] = np.max(np.abs(sol.y[0]))

    # === Affichage du diagramme ===
    plt.figure(figsize=(10, 6))
    plt.contourf(omegas, As, amplitudes_max, levels=30, cmap="plasma")
    plt.colorbar(label="Amplitude maximale de $\\theta_1$")
    plt.xlabel("Fréquence d'activation $\\omega$")
    plt.ylabel("Amplitude d'activation $A$")
    plt.title("Diagramme de stabilité avec inhibition lissée pour N = {}".format(N))
    plt.show()