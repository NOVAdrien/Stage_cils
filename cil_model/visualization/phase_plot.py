# visualization/phase_plot.py

import matplotlib.pyplot as plt
import numpy as np

def plot_phase(solution, N, mode="single", ncols=3):
    """
    Plot des portraits de phase (θ_i vs θ̇_i).

    :param solution: objet solution de solve_ivp
    :param N: nombre de barres
    :param mode: "single" (toutes les courbes sur 1 graphe) ou "subplots" (1 graphe par θ_i)
    :param ncols: nombre de colonnes pour les subplots
    """
    angles = solution.y[:N]         # θ_i(t)
    speeds = solution.y[N:2*N]      # θ̇_i(t)

    if mode == "single":
        plt.figure(figsize=(8, 6))
        for i in range(N):
            plt.plot(angles[i], speeds[i], label=f'θ_{i+1} vs θ̇_{i+1}')
        plt.xlabel('theta_i (rad)')
        plt.ylabel('theta_dot_i(rad/s)')
        plt.title(f'Portraits de phase (N = {N})')
        plt.legend()
        plt.grid()
        plt.show()

    elif mode == "subplots":
        nrows = (N + ncols - 1) // ncols  # arrondi vers le haut
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3), sharex=False, sharey=False)

        # Si axes est 1D → le rendre 2D pour simplifier la boucle
        axes = np.array(axes).reshape(nrows, ncols)

        for idx in range(N):
            row, col = divmod(idx, ncols)
            ax = axes[row, col]
            ax.plot(angles[idx], speeds[idx], label=f'θ_{idx+1} vs θ̇_{idx+1}')
            ax.set_xlabel('θ_i (rad)')
            ax.set_ylabel('θ̇_i (rad/s)')
            ax.grid()
            ax.legend()

        # Supprimer les subplots vides
        for idx in range(N, nrows * ncols):
            row, col = divmod(idx, ncols)
            fig.delaxes(axes[row, col])

        fig.suptitle(f'Portraits de phase (N = {N})', fontsize=16)
        plt.tight_layout()
        plt.show()

    else:
        raise ValueError(f"Mode inconnu: {mode}. Utiliser 'single' ou 'subplots'.")
