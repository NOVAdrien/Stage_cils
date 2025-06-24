# visualization/animate.py

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def animate_chain(solution, N, lengths, save=False, filename="animation.mp4"):
    """
    Animation du mouvement de la chaîne de barres rigides (cil).

    :param solution: objet solution de solve_ivp
    :param N: nombre de barres
    :param lengths: liste des longueurs des barres
    :param save: bool, si True → sauvegarde en .mp4
    :param filename: nom du fichier de sortie .mp4
    """
    t = solution.t
    angles = solution.y[:N]  # θ_i(t)

    # === Calcul des positions cumulées ===
    # Pour chaque frame (chaque t[k]), on calcule (x0,y0), (x1,y1), ..., (xN,yN)
    x_data = []
    y_data = []

    for k in range(len(t)):
        x_positions = [0]
        y_positions = [0]
        x_current, y_current = 0, 0
        cumulative_angle = 0

        for i in range(N):
            cumulative_angle += angles[i][k]
            dx = lengths[i] * np.sin(cumulative_angle)
            dy = lengths[i] * np.cos(cumulative_angle)

            x_current += dx
            y_current += dy

            x_positions.append(x_current)
            y_positions.append(y_current)

        x_data.append(x_positions)
        y_data.append(y_positions)

    # === Setup du plot ===
    fig, ax = plt.subplots(figsize=(6, 6))
    line, = ax.plot([], [], 'o-', lw=3)

    # Fixer les limites en fonction de la longueur totale
    L_totale = sum(lengths)
    ax.set_xlim(-L_totale * 1.2, L_totale * 1.2)
    ax.set_ylim(0, L_totale * 1.5)
    ax.set_aspect('equal')
    ax.grid()
    ax.set_title(f'Animation du cil (N = {N})')

    # === Fonction d'update pour l'animation ===
    def update(frame):
        line.set_data(x_data[frame], y_data[frame])
        return line,

    ani = animation.FuncAnimation(fig, update, frames=len(t), interval=0.002, blit=True)

    # === Sauvegarde si demandé ===
    if save:
        ani.save(filename, fps=30, dpi=200)
        print(f"Animation sauvegardée sous : {filename}")

    plt.show()
