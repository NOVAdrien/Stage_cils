# visualization/plot_angles.py

import matplotlib.pyplot as plt
import numpy as np

def plot_angles(solution, N, mode="single", ncols=3):
    """
    Plot des angles θ_i(t).

    :param solution: objet solution de solve_ivp
    :param N: nombre de barres
    :param mode: "single" (toutes les courbes sur 1 graphe) ou "subplots" (1 graphe par θ_i)
    :param ncols: nombre de colonnes pour les subplots
    """
    t = solution.t
    angles = solution.y[:N]  # les N premières lignes de y = les θ_i(t)

    if mode == "single":
        plt.figure(figsize=(8, 6))
        for i in range(N):
            plt.plot(t, angles[i], label=f'θ_{i+1}(t)')
        plt.xlabel('Temps (s)')
        plt.ylabel('Angle (rad)')
        plt.title(f'Évolution des angles (N = {N})')
        plt.legend()
        plt.grid()
        plt.show()

    elif mode == "subplots":
        nrows = (N + ncols - 1) // ncols  # arrondi vers le haut
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3), sharex=True)

        # Si axes est 1D → le rendre 2D pour simplifier la boucle
        axes = np.array(axes).reshape(nrows, ncols)

        for idx in range(N):
            row, col = divmod(idx, ncols)
            ax = axes[row, col]
            ax.plot(t, angles[idx], label=f'θ_{idx+1}(t)')
            ax.set_xlabel('Temps (s)')
            ax.set_ylabel('Angle (rad)')
            ax.grid()
            ax.legend()

        # Supprimer les subplots vides
        for idx in range(N, nrows * ncols):
            row, col = divmod(idx, ncols)
            fig.delaxes(axes[row, col])

        fig.suptitle(f'Évolution des angles (N = {N})', fontsize=16)
        plt.tight_layout()
        plt.show()

    else:
        raise ValueError(f"Mode inconnu: {mode}. Utiliser 'single' ou 'subplots'.")
