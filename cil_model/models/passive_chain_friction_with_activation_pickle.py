# models/passive_chain_friction_with_activation_pickle.py

import numpy as np
import sympy as sp
from sympy import diff
from sympy.physics.mechanics import dynamicsymbols
from lagrangian.energy_builder import (
    energie_cinetique_translation,
    energie_cinetique_rotation,
    energie_potentielle_ressorts_activation
)
from lagrangian.dissipation_builder import (
    fonction_dissipation,
    calcul_Q_non_cons
)
from lagrangian.frottements_solver import (
    compute_M_list,
    ajoute_frottements_au_rhs
)
from lagrangian.lagrangian_solver import generate_lagrangian_rhs
from scipy.integrate import solve_ivp

import time
import os
import dill as pickle

from .base_system import BaseSystem

class PassiveChainFrictionWithActivationPickle(BaseSystem):
    """
    Système de N barres rigides AVEC frottements et activation interne.
    """

    def __init__(self, N, masses, lengths, widths, stiffnesses, amplitudes, omega, k_spatial, phi, alpha, gamma, initial_conditions):
        
        start_time_symb = time.time()
        super().__init__(N)
        self.masses = masses
        self.lengths = lengths
        self.widths = widths
        self.stiffnesses = stiffnesses
        self.initial_conditions = initial_conditions

        # Symboles de frottements
        self.alpha = sp.symbols(f'alpha1:{self.N+1}')
        self.gamma = sp.symbols(f'gamma1:{self.N+1}')
        self.alpha_values = alpha
        self.gamma_values = gamma

        # Symboles de activations
        self.a = sp.symbols(f'amlpitude1:{self.N+1}')
        self.o = sp.symbols(f'omega1:{self.N+1}')
        self.k_sp = sp.symbols(f'k_spatial1:{self.N+1}')
        self.phi_symb = sp.symbols(f'phi1:{self.N+1}')

        # Stockage des valeurs numériques
        self.amplitudes = amplitudes
        self.omega = omega
        self.k_spatial = k_spatial
        self.phi = phi

        # Variables symboliques
        self.t = sp.symbols('t')
        self.theta = dynamicsymbols(f'theta1:{N+1}')
        self.theta_dot = [diff(th, self.t) for th in self.theta]

        # Symboles physiques
        self.m = sp.symbols(f'm1:{N+1}')
        self.l = sp.symbols(f'l1:{N+1}')
        self.r = sp.symbols(f'r1:{N+1}')
        self.k = sp.symbols(f'k1:{N+1}')
        print("Initialisation symbolique : ", time.time() - start_time_symb)

        start_time = time.time()
        # On reconstruit le Lagrangien
        self.L = self._build_lagrangian()
        print("build_lagrangian : ", time.time() - start_time)

        # === Fichier de cache pour rhs_modifie déjà numérisé ===
        cache_file = f"rhs_modifie_num_N{N}.pkl"
        start_time_ifelse = time.time()


        if os.path.exists(cache_file) and os.path.getsize(cache_file) > 0:
            print("Lecture ...")
            # === Chargement direct du rhs_modifie avec subs déjà faits ===
            with open(cache_file, "rb") as f:
                rhs_with_friction_and_activation_pickle = pickle.load(f)
            print(f"rhs_modifie (numérique) chargé depuis {cache_file}")
            print("cache :", time.time() - start_time_ifelse)

            # args = [self.t] + list(self.theta) + list(self.theta_dot) \
            #     + list(self.m) + list(self.l) + list(self.r) + list(self.k) + list(self.a) + list(self.o) + list(self.k_sp) + list(self.phi_symb) \
            #     + list(self.alpha) + list(self.gamma)
            
            # self._rhs_func = sp.lambdify(args, rhs_with_friction_and_activation_pickle)

        else:
            print("Calcul symbolique + subs...")

            # # On reconstruit le Lagrangien
            # self.L = self._build_lagrangian()

            # Construire la dissipation
            self.compute_dissipation()
            self.compute_Q()
            self.compute_M()

            # Construire le RHS avec frottements et activation
            self.LM, self._rhs_func = generate_lagrangian_rhs(
                self.L, self.theta, self.t, [self.m, self.l, self.r, self.k, self.a, self.o, self.k_sp, self.phi_symb]
            )

            # On calcule le RHS symbolique
            rhs_with_friction_and_activation_pickle = self.LM.rhs()

            # On ajoute les frottements :
            rhs_with_friction_and_activation_pickle = ajoute_frottements_au_rhs(rhs_with_friction_and_activation_pickle, self.Q_non_cons, self.M_list, self.N)

            # Sauvegarde du rhs_modifie numérisé
            with open(cache_file, "wb") as f:
                pickle.dump(rhs_with_friction_and_activation_pickle, f)
            print(f"rhs_modifie numérisé sauvegardé dans {cache_file}")
            print("Calcul else :", time.time() - start_time_ifelse)

        # Lambdify final :
        args = [self.t] + list(self.theta) + list(self.theta_dot) \
            + list(self.m) + list(self.l) + list(self.r) + list(self.k) + list(self.a) + list(self.o) + list(self.k_sp) + list(self.phi_symb) \
            + list(self.alpha) + list(self.gamma)

        self._rhs_func = sp.lambdify(args, rhs_with_friction_and_activation_pickle)
            
    def _build_lagrangian(self):
        self.Ec_trans = energie_cinetique_translation(
            self.N, self.theta, self.t, self.m, self.l
        )
        self.Ec_rot = energie_cinetique_rotation(
            self.N, self.theta_dot, self.m, self.l, self.r
        )
        self.Ec = self.Ec_trans + self.Ec_rot

        self.Ep = energie_potentielle_ressorts_activation(
            self.N, self.theta, self.t, self.k, self.l, self.amplitudes, self.omega, self.k_spatial, self.phi
        )

        L = self.Ec - self.Ep
        return L
    
    def compute_dissipation(self):
        self.D = fonction_dissipation(
            self.N, self.theta, self.theta_dot, self.t,
            self.alpha, self.gamma, self.m, self.l
        )

    def compute_Q(self):
        self.Q_non_cons = calcul_Q_non_cons(self.N, self.D, self.theta_dot, self.t)

    def compute_M(self):
        self.M_list = compute_M_list(self.N, self.L, self.theta_dot)


    def solve(self, t_span, t_eval=None, method='RK45'):
        """
        Résout les équations du système AVEC frottements et activation interne.
        """
        def rhs(t_num, y):
            theta_vals = y[:self.N]
            theta_dot_vals = y[self.N:]
            
            # Construction des inputs
            inputs = [t_num] + list(theta_vals) + list(theta_dot_vals) \
                + list(self.masses) + list(self.lengths) + list(self.widths) + list(self.stiffnesses) + list(self.amplitudes) + list(self.omega) + list(self.k_spatial) + list(self.phi) \
                + list(self.alpha_values) + list(self.gamma_values)
            
            out = self._rhs_func(*inputs)
            return np.array(out, dtype=float).flatten()
        
        angles_init, speeds_init = self.initial_conditions
        print("angles_init = ", angles_init)
        print("speeds_init = ", speeds_init)
        y0 = np.concatenate([angles_init, speeds_init])
        
        # Résolution avec solve_ivp
        sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, method=method, rtol=1e-9, atol=1e-9)

        return sol