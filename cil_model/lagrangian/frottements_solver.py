 # lagrangian/frottements_solver.py

import sympy as sp
from sympy import diff

def compute_M_list(N, L, theta_dot):
    """
    Calcule les M_j = ∂²L / ∂theta_dot[j]² pour chaque j.
    """
    M_list = []
    for j in range(N):
        Mj = diff(diff(L, theta_dot[j]), theta_dot[j])
        M_list.append(Mj)
    return M_list

def ajoute_frottements_au_rhs(rhs, Q_non_cons, M_list, N):
    """
    Modifie le RHS en ajoutant les termes de frottement :
    rhs[N + j] += Q_j / M_j
    """
    rhs_modifie = list(rhs)  # On fait une copie pour ne pas modifier l'original
    for j in range(N):
        rhs_modifie[N + j] += Q_non_cons[j] / M_list[j]
    return rhs_modifie
