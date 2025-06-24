# models/passive_chain.py

import numpy as np
import sympy as sp
from sympy import diff
from sympy.physics.mechanics import dynamicsymbols
from scipy.integrate import solve_ivp
from lagrangian.energy_builder import (
    energie_cinetique_translation,
    energie_cinetique_rotation,
    energie_potentielle_ressorts
)
from lagrangian.lagrangian_solver import generate_lagrangian_rhs



from .base_system import BaseSystem

class PassiveChain(BaseSystem):
    """
    Système de N barres rigides passives (ressorts sur les pivots).
    """

    def __init__(self, N, masses, lengths, widths, stiffnesses, initial_conditions):
        """
        Initialise le système passif.

        :param N: Nombre de barres
        :param masses: Liste des masses [m1, m2, ..., mN]
        :param lengths: Liste des longueurs [l1, l2, ..., lN]
        :param widths: Liste des demi-largeurs [r1, r2, ..., rN]
        :param stiffnesses: Liste des raideurs de ressort [k1, k2, ..., kN]
        :param initial_conditions: Tuple (angles_init, vitesses_init)
        """
        super().__init__(N)
        self.masses = masses
        self.lengths = lengths
        self.widths = widths
        self.stiffnesses = stiffnesses
        self.initial_conditions = initial_conditions

        # Variables symboliques
        self.t = sp.symbols('t')
        self.theta = dynamicsymbols(f'theta1:{N+1}')
        self.theta_dot = [diff(th, self.t) for th in self.theta]

        # Symboles physiques
        self.m = sp.symbols(f'm1:{N+1}')
        self.l = sp.symbols(f'l1:{N+1}')
        self.r = sp.symbols(f'r1:{N+1}')
        self.k = sp.symbols(f'k1:{N+1}')

        # Construction du Lagrangien
        self.L = self._build_lagrangian()

        self.LM, self._rhs_func = generate_lagrangian_rhs(
            self.L, self.theta, self.t, [self.m, self.l, self.r, self.k]
        )


    def _build_lagrangian(self):
        self.Ec_trans = energie_cinetique_translation(
            self.N, self.theta, self.t, self.m, self.l
        )
        self.Ec_rot = energie_cinetique_rotation(
            self.N, self.theta_dot, self.m, self.l, self.r
        )
        self.Ec = self.Ec_trans + self.Ec_rot

        self.Ep = energie_potentielle_ressorts(
            self.N, self.theta, self.k
        )

        L = self.Ec - self.Ep
        return L



    def solve(self, t_span, t_eval=None, method='RK45'):
        """
        Résout les équations du système.

        :param t_span: Intervalle de temps (t0, tf)
        :param t_eval: Points de temps pour l'échantillonnage
        :param method: Méthode d'intégration
        :return: solution de solve_ivp
        """
        def rhs(t_num, y):
            theta_vals = y[:self.N]
            theta_dot_vals = y[self.N:]

            # Construction de l'entrée de la fonction RHS
            inputs = [t_num] + list(theta_vals) + list(theta_dot_vals) \
                     + list(self.masses) + list(self.lengths) + list(self.widths) + list(self.stiffnesses)

            out = self._rhs_func(*inputs)
            return np.array(out, dtype=float).flatten()

        angles_init, speeds_init = self.initial_conditions
        y0 = np.concatenate([angles_init, speeds_init])

        sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, method=method, rtol=1e-9, atol=1e-9)
        return sol
