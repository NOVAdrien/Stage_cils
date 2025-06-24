# lagrangian/energy_builder.py

import sympy as sp
from sympy import sin, cos, tanh, diff

def energie_cinetique_translation(N, theta, t, m, l):
    """
    Energie cinétique de translation des centres de masse.
    """
    x, y = [], []
    for i in range(N):
        if i == 0:
            xi = -(l[i] / 2) * sin(theta[i])
            yi =  (l[i] / 2) * cos(theta[i])
        else:
            xi = x[i-1] - (l[i-1]/2) * sin(theta[i-1]) - (l[i]/2) * sin(theta[i])
            yi = y[i-1] + (l[i-1]/2) * cos(theta[i-1]) + (l[i]/2) * cos(theta[i])
        x.append(xi)
        y.append(yi)

    Ec_trans = 0
    for i in range(N):
        vx = diff(x[i], t)
        vy = diff(y[i], t)
        Ec_trans += (1/2) * m[i] * (vx**2 + vy**2)
    return Ec_trans

def energie_cinetique_rotation(N, theta_dot, m, l, r):
    """
    Energie cinétique de rotation.
    """
    Ec_rot = 0
    for i in range(N):
        I_i = m[i] * (l[i]**2 / 12 + r[i]**2 / 3)
        Ec_rot += (1/2) * I_i * theta_dot[i]**2
    return Ec_rot

def energie_potentielle_ressorts(N, theta, k):
    """
    Energie potentielle des ressorts sur les pivots.
    """
    Ep = 0
    for i in range(N):
        theta_prev = 0 if i == 0 else theta[i - 1]
        Ep += (1/2) * k[i] * (theta[i] - theta_prev)**2
    return Ep

def energie_potentielle_ressorts_activation(N, theta, t, k, l, A, omega, k_spatial, phi):
    """
    Energie potentielle avec activation interne des cils :
    Ep_i = 1/2 * k_i * (theta_i - theta_{i-1} - theta_activ_i(t))^2
    """
    Ep = 0
    s_i = 0

    for i in range(N):
        s_i += l[i]
        if i == 0:
            theta_diff = theta[i]
        else:
            theta_diff = theta[i] - theta[i-1]

        theta_activ_i = A[i] * sin(omega[i] * t - k_spatial[i] * s_i + phi[i])

        Ep += (1/2) * k[i] * (theta_diff - theta_activ_i)**2

    return Ep

def energie_potentielle_ressorts_activation_critique(N, theta, t, k, l, A, omega, k_spatial, phi, eps, tc):
    """
    Energie potentielle avec activation interne des cils :
    Ep_i = 1/2 * k_i * (theta_i - theta_{i-1} - theta_activ_i(t))^2
    """
    Ep = 0
    s_i = 0
    print("Theta = ", theta)
    print("Theta_crit = ", tc)
    print("Epsilon = ", eps)

    for i in range(N):
        s_i += l[i]
        if i == 0:
            theta_diff = theta[i]
        else:
            theta_diff = theta[i] - theta[i-1]

        inhibition_i = 0.5*(1 + tanh(-(theta[i] - tc[i])**2 / eps[i]))
        print(f"Inhibition pour cil {i}: {inhibition_i}")
        theta_activ_i = inhibition_i * A[i] * sin(omega[i] * t - k_spatial[i] * s_i + phi[i])

        Ep += (1/2) * k[i] * (theta_diff - theta_activ_i)**2

    return Ep