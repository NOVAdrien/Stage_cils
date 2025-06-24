# models/passive_chain_friction.py

import numpy as np
import sympy as sp

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

from .passive_chain import PassiveChain

class PassiveChainWithFriction(PassiveChain):
    """
    Système de N barres rigides AVEC frottements.
    """

    def __init__(self, N, masses, lengths, widths, stiffnesses, alpha, gamma, initial_conditions):
        super().__init__(N, masses, lengths, widths, stiffnesses, initial_conditions)
        
        # Symboles de frottements
        self.alpha = sp.symbols(f'alpha1:{self.N+1}')
        self.gamma = sp.symbols(f'gamma1:{self.N+1}')
        
        # Stockage des valeurs numériques
        self.alpha_values = alpha
        self.gamma_values = gamma
        
        # Construire la dissipation
        self.compute_dissipation()
        self.compute_Q()
        self.compute_M()
        
        # Construire le RHS avec frottements
        self._build_rhs_with_friction()
    
    def compute_dissipation(self):
        self.D = fonction_dissipation(
            self.N, self.theta, self.theta_dot, self.t,
            self.alpha, self.gamma, self.m, self.l
        )
    
    def compute_Q(self):
        self.Q_non_cons = calcul_Q_non_cons(self.N, self.D, self.theta_dot, self.t)
    
    def compute_M(self):
        self.M_list = compute_M_list(self.N, self.L, self.theta_dot)

    
    def _build_rhs_with_friction(self):
        # On repart de la version SANS frottement :
        self.LM, rhs_no_friction_func = generate_lagrangian_rhs(
            self.L, self.theta, self.t, [self.m, self.l, self.r, self.k]
        )
        
        # On récupère le RHS symbolique
        rhs_no_friction = self.LM.rhs()
        
        # On ajoute les frottements :
        rhs_with_friction = ajoute_frottements_au_rhs(rhs_no_friction, self.Q_non_cons, self.M_list, self.N)
        
        # Lambdify final :
        args = [self.t] + list(self.theta) + list(self.theta_dot) \
             + list(self.m) + list(self.l) + list(self.r) + list(self.k) \
             + list(self.alpha) + list(self.gamma)
        
        self._rhs_func = sp.lambdify(args, rhs_with_friction)
    
    def solve(self, t_span, t_eval=None, method='RK45'):
        """
        Résout les équations du système AVEC frottements.
        """
        def rhs(t_num, y):
            theta_vals = y[:self.N]
            theta_dot_vals = y[self.N:]
            
            # Construction des inputs
            inputs = [t_num] + list(theta_vals) + list(theta_dot_vals) \
                   + list(self.masses) + list(self.lengths) + list(self.widths) + list(self.stiffnesses) \
                   + list(self.alpha_values) + list(self.gamma_values)
            
            out = self._rhs_func(*inputs)
            return np.array(out, dtype=float).flatten()
        
        angles_init, speeds_init = self.initial_conditions
        y0 = np.concatenate([angles_init, speeds_init])
        
        sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, method=method, rtol=1e-9, atol=1e-9)
        
        return sol
