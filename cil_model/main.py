# main_test_all.py

from models.passive_chain import PassiveChain
from models.passive_chain_friction import PassiveChainWithFriction
from models.passive_chain_friction_with_activation import PassiveChainFrictionWithActivation
from models.passive_chain_friction_with_activation_pickle import PassiveChainFrictionWithActivationPickle
from models.passive_chain_friction_with_activation_pickle_crit import PassiveChainFrictionWithActivationPickleCrit
from visualization.plot_angles import plot_angles
from visualization.phase_plot import plot_phase
from visualization.animate import animate_chain
from visualization.energy_plot import plot_energy, plot_energy_crit
from visualization.spectrum_plot import plot_spectrum
from visualization.plot_mathieu import plot_mathieu

import sympy as sp
import numpy as np
import time

# === Settings === #
start_time = time.time()
N = 2

masses = [1.0] * N
lengths = [1.0] * N
widths  = [0.2] * N
stiffnesses = [10.0] * N
amplitudes = [0.1] * N
omega = [2.0] * N
k_spatial = [1.0] * N
phi = [0.0] * N

# === Friction === #
alpha = [0.1] * N
gamma = [0.1] * N

# === Initial conditions === #
angles_init = [0.5] * N
speeds_init = [0.0] * N
initial_conditions = (angles_init, speeds_init)

# === Time span and evaluation points === #
t_span = (0, 20)
t_eval = np.linspace(0, 20, 10000)

# Inhibition
epsilon = [0.1] * N
theta_crit = [0.6] * N

print("Initializing parameters: ", time.time() - start_time)










# === MODELE 1 : PassiveChain (SANS frottement) ===
start_time = time.time()
print("=== PassiveChain (sans frottement) ===")

chain = PassiveChain(
    N, masses, lengths, widths, stiffnesses,
    initial_conditions
)

sol = chain.solve(t_span, t_eval=t_eval, method='RK45')

# Lambdify des énergies
Ec_func = sp.lambdify(
    [chain.t] + list(chain.theta) + list(chain.theta_dot) \
    + list(chain.m) + list(chain.l) + list(chain.r) + list(chain.k),
    chain.Ec
)

Ep_func = sp.lambdify(
    [chain.t] + list(chain.theta) + list(chain.theta_dot) \
    + list(chain.m) + list(chain.l) + list(chain.r) + list(chain.k),
    chain.Ep
)

print("Total time taken for PassiveChain: ", time.time() - start_time)

# Plots
plot_angles(sol, N, mode="single")
plot_phase(sol, N, mode="single")
plot_energy(sol, Ec_func, Ep_func, t_eval, N, chain, title="Énergie mécanique - PassiveChain")
plot_spectrum(sol, N, t_eval)
animate_chain(sol, N, lengths, save=False)










# === MODELE 2 : PassiveChainWithFriction ===
start_time = time.time()
print("\n=== PassiveChainWithFriction ===")

chain_friction = PassiveChainWithFriction(
    N, masses, lengths, widths, stiffnesses,
    alpha, gamma,
    initial_conditions
)

sol_friction = chain_friction.solve(t_span, t_eval=t_eval, method='RK45')

# Lambdify des énergies
Ec_func_friction = sp.lambdify(
    [chain_friction.t] + list(chain_friction.theta) + list(chain_friction.theta_dot) \
    + list(chain_friction.m) + list(chain_friction.l) + list(chain_friction.r) + list(chain_friction.k),
    chain_friction.Ec
)

Ep_func_friction = sp.lambdify(
    [chain_friction.t] + list(chain_friction.theta) + list(chain_friction.theta_dot) \
    + list(chain_friction.m) + list(chain_friction.l) + list(chain_friction.r) + list(chain_friction.k),
    chain_friction.Ep
)

print("Total time taken for PassiveChainWithFriction: ", time.time() - start_time)

# Plots
plot_angles(sol_friction, N, mode="single")
plot_phase(sol_friction, N, mode="single")
plot_energy(sol_friction, Ec_func_friction, Ep_func_friction, t_eval, N, chain_friction, title="Énergie mécanique - PassiveChainWithFriction")
plot_spectrum(sol_friction, N, t_eval)
animate_chain(sol_friction, N, lengths, save=False)










# === MODELE 3 : PassiveChainWithFrictionAndActivations ===
start_time = time.time()
print("\n=== PassiveChainWithFrictionAndActivation ===")

chain_activation = PassiveChainFrictionWithActivation(
    N, masses, lengths, widths, stiffnesses, amplitudes, omega, k_spatial, phi,
    alpha, gamma,
    initial_conditions
)

sol_activation = chain_activation.solve(t_span, t_eval=t_eval, method='RK45')

# Lambdify des énergies
Ec_func_activation = sp.lambdify(
    [chain_activation.t] + list(chain_activation.theta) + list(chain_activation.theta_dot) \
    + list(chain_activation.m) + list(chain_activation.l) + list(chain_activation.r) + list(chain_activation.k),
    chain_activation.Ec
)

Ep_func_activation = sp.lambdify(
    [chain_activation.t] + list(chain_activation.theta) + list(chain_activation.theta_dot) \
    + list(chain_activation.m) + list(chain_activation.l) + list(chain_activation.r) + list(chain_activation.k),
    chain_activation.Ep
)

print("Total time taken for PassiveChainWithFrictionAndActivation: ", time.time() - start_time)

# Plots
plot_angles(sol_activation, N, mode="single")
plot_phase(sol_activation, N, mode="single")
plot_energy(sol_activation, Ec_func_activation, Ep_func_activation, t_eval, N, chain_activation, title="Énergie mécanique - PassiveChainWithFrictionAndActivation")
plot_spectrum(sol_activation, N, t_eval)
animate_chain(sol_activation, N, lengths, save=False)










# === MODELE 4 : Test pickle for PassiveChainWithFrictionAndActivationsPickle ===
start_time = time.time()
print("\n=== Test pickle for PassiveChainWithFrictionAndActivationPickle ===")

chain_activation = PassiveChainFrictionWithActivationPickle(
    N, masses, lengths, widths, stiffnesses, amplitudes, omega, k_spatial, phi,
    alpha, gamma,
    initial_conditions
)

sol_activation = chain_activation.solve(t_span, t_eval=t_eval, method='RK45')

# Lambdify des énergies
Ec_func_activation = sp.lambdify(
    [chain_activation.t] + list(chain_activation.theta) + list(chain_activation.theta_dot) \
    + list(chain_activation.m) + list(chain_activation.l) + list(chain_activation.r) + list(chain_activation.k),
    chain_activation.Ec
)

Ep_func_activation = sp.lambdify(
    [chain_activation.t] + list(chain_activation.theta) + list(chain_activation.theta_dot) \
    + list(chain_activation.m) + list(chain_activation.l) + list(chain_activation.r) + list(chain_activation.k),
    chain_activation.Ep
)

print("Total time taken for PassiveChainWithFrictionAndActivationPickle: ", time.time() - start_time)

# Plots
plot_angles(sol_activation, N, mode="single")
plot_phase(sol_activation, N, mode="single")
plot_energy(sol_activation, Ec_func_activation, Ep_func_activation, t_eval, N, chain_activation, title="Énergie mécanique - PassiveChainWithFrictionAndActivationPickle")
plot_spectrum(sol_activation, N, t_eval)
animate_chain(sol_activation, N, lengths, save=False)










# === MODELE 5 : Test pickle for PassiveChainWithFrictionAndActivationsPickleCrit ===
start_time = time.time()
print("\n=== Test pickle for PassiveChainWithFrictionAndActivationPickleCrit ===")

chain_crit = PassiveChainFrictionWithActivationPickleCrit(
    N, masses, lengths, widths, stiffnesses, amplitudes, omega, k_spatial, phi,
    alpha, gamma,
    initial_conditions,
    epsilon, theta_crit
)

sol_crit = chain_crit.solve(t_span, t_eval=t_eval, method='RK45')

# Lambdify des énergies
Ec_func_crit = sp.lambdify(
    [chain_crit.t] + list(chain_crit.theta) + list(chain_crit.theta_dot) \
    + list(chain_crit.m) + list(chain_crit.l) + list(chain_crit.r) + list(chain_crit.k) + list(chain_crit.eps) + list(chain_crit.tc),
    chain_crit.Ec
)

Ep_func_crit = sp.lambdify(
    [chain_crit.t] + list(chain_crit.theta) + list(chain_crit.theta_dot) \
    + list(chain_crit.m) + list(chain_crit.l) + list(chain_crit.r) + list(chain_crit.k) + list(chain_crit.eps) + list(chain_crit.tc),
    chain_crit.Ep
)

print("Total time taken for PassiveChainWithFrictionAndActivationPickleCrit: ", time.time() - start_time)

# Plots
plot_angles(sol_crit, N, mode="single")
plot_phase(sol_crit, N, mode="single")
plot_energy_crit(sol_crit, Ec_func_crit, Ep_func_crit, t_eval, N, chain_crit, title="Énergie mécanique - PassiveChainWithFrictionAndActivationPickleCrit")
plot_spectrum(sol_crit, N, t_eval)
animate_chain(sol_crit, N, lengths, save=False)










# === MODELE 6 (BONUS) : Mathieu diagram ===
plot_mathieu(masses, lengths, widths, stiffnesses, k_spatial, phi, alpha, gamma, initial_conditions, t_span, N, t_eval)