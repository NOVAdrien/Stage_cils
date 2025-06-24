# lagrangian/dissipation_builder.py

import sympy as sp
from sympy import diff, sin, cos

def fonction_dissipation(N, theta, theta_dot, t, alpha, gamma, m, l):
    """
    Fonction de dissipation D :
    D = sum_i (1/2) * alpha_i * v_G_i^2 + sum_i (1/2) * gamma_i * theta_dot_i^2
    """
    # Calcul des positions des centres de masse
    x, y = [], []
    for i in range(N):
        if i == 0:
            xi = -(l[i] / 2) * sp.sin(theta[i])
            yi =  (l[i] / 2) * sp.cos(theta[i])
        else:
            xi = x[i-1] - (l[i-1]/2) * sp.sin(theta[i-1]) - (l[i]/2) * sp.sin(theta[i])
            yi = y[i-1] + (l[i-1]/2) * sp.cos(theta[i-1]) + (l[i]/2) * sp.cos(theta[i])
        x.append(xi)
        y.append(yi)

    # Dissipation translation
    D_trans = 0
    for i in range(N):
        vx = diff(x[i], t)
        vy = diff(y[i], t)
        v_squared = vx**2 + vy**2
        D_trans += (1/2) * alpha[i] * v_squared

    # Dissipation rotation
    D_rot = 0
    for i in range(N):
        D_rot += (1/2) * gamma[i] * theta_dot[i]**2

    # Total
    D_total = D_trans + D_rot
    return D_total

def calcul_Q_non_cons(N, D, theta_dot, t):
    """
    Calcule les forces généralisées de dissipation :
    Q_j = - ∂D / ∂theta_dot_j
    """
    Q_non_cons = []
    for j in range(N):
        Qj = - diff(D, theta_dot[j])
        Q_non_cons.append(Qj)
    return Q_non_cons
