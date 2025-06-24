# lagrangian/lagrangian_solver.py

import sympy as sp
from sympy.physics.mechanics import LagrangesMethod

def generate_lagrangian_rhs(L, theta, t, symboles_parametres):
    """
    Génère le LagrangesMethod et la fonction RHS lambdifiée.

    :param L: Lagrangien (sympy expression)
    :param theta: liste des variables theta_i(t)
    :param t: symbole du temps
    :param symboles_parametres: liste de listes de symboles physiques [m, l, r, k, ...]
    :return: (LagrangesMethod, rhs_func)
    """
    LM = LagrangesMethod(L, theta)
    LM.form_lagranges_equations()

    # Construire la liste des arguments pour le lambdify
    args = [t] + list(theta) + [th.diff(t) for th in theta]
    for symbole_list in symboles_parametres:
        args += list(symbole_list)

    # Lambdification
    rhs_func = sp.lambdify(args, LM.rhs())

    return LM, rhs_func
