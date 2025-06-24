# visualization/energy_plot.py

import numpy as np
import matplotlib.pyplot as plt

def plot_energy(sol, Ec_func, Ep_func, t_eval, N, chain, title="Énergie mécanique"):
    """
    Trace l'énergie mécanique totale E_c + E_p en fonction du temps.
    """
    Ec_vals = []
    Ep_vals = []
    for i, t_num in enumerate(t_eval):
        theta_vals = sol.y[:N, i]
        theta_dot_vals = sol.y[N:, i]
        
        inputs = [t_num] + list(theta_vals) + list(theta_dot_vals) \
               + list(chain.masses) + list(chain.lengths) + list(chain.widths) + list(chain.stiffnesses)
        
        Ec = Ec_func(*inputs)
        Ep = Ep_func(*inputs)
        
        Ec_vals.append(Ec)
        Ep_vals.append(Ep)
    
    Ec_vals = np.array(Ec_vals)
    Ep_vals = np.array(Ep_vals)
    E_mec = Ec_vals + Ep_vals
    
    plt.figure(figsize=(8, 5))
    plt.plot(t_eval, E_mec, label=r"$E_{\mathrm{mécanique}}$")
    plt.plot(t_eval, Ec_vals, '--', label=r"$E_c$")
    plt.plot(t_eval, Ep_vals, '--', label=r"$E_p$")
    plt.xlabel("Temps t")
    plt.ylabel("Énergie")
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.show()

def plot_energy_crit(sol, Ec_func, Ep_func, t_eval, N, chain, title="Énergie mécanique"):
    """
    Trace l'énergie mécanique totale E_c + E_p en fonction du temps.
    """
    Ec_vals = []
    Ep_vals = []
    for i, t_num in enumerate(t_eval):
        theta_vals = sol.y[:N, i]
        theta_dot_vals = sol.y[N:, i]
        
        inputs = [t_num] + list(theta_vals) + list(theta_dot_vals) \
               + list(chain.masses) + list(chain.lengths) + list(chain.widths) + list(chain.stiffnesses) + list(chain.epsilon) + list(chain.theta_crit)
        
        Ec = Ec_func(*inputs)
        Ep = Ep_func(*inputs)
        
        Ec_vals.append(Ec)
        Ep_vals.append(Ep)
    
    Ec_vals = np.array(Ec_vals)
    Ep_vals = np.array(Ep_vals)
    E_mec = Ec_vals + Ep_vals
    
    plt.figure(figsize=(8, 5))
    plt.plot(t_eval, E_mec, label=r"$E_{\mathrm{mécanique}}$")
    plt.plot(t_eval, Ec_vals, '--', label=r"$E_c$")
    plt.plot(t_eval, Ep_vals, '--', label=r"$E_p$")
    plt.xlabel("Temps t")
    plt.ylabel("Énergie")
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.show()