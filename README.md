# Stage_cils

cil_model/
│
├── main.py                            # Script principal (lance les simulations, affiche)
├── config.py                          # Paramètres globaux et presets de systèmes
│
├── models/
│   ├── __init__.py
│   ├── base_system.py                 # Classe abstraite : structure commune
│   ├── passive_chain.py              # Système de N barres sans moteur ni frottement
│   ├── active_chain.py               # Avec moteurs (structures actives)
│   └── friction_mixins.py            # Mixins ou classes utilitaires pour les frottements
│
├── lagrangian/
│   ├── __init__.py
│   ├── energy_builder.py             # Calcul des énergies : Ec_trans, Ec_rot, Ep
│   ├── lagrangian_solver.py          # Génère et résout les équations d’Euler-Lagrange
│   └── symbolic_helpers.py           # Fonctions sympy réutilisables
│
├── simulation/
│   ├── __init__.py
│   ├── integrator.py                 # Appel à solve_ivp, configuration du système
│   └── analyser.py                   # Portraits de phase, spectre, stabilité...
│
├── visualization/
│   ├── __init__.py
│   ├── plot_angles.py                # Traces θ₁(t), θ₂(t), etc.
│   ├── phase_plot.py                 # Portraits de phase
│   └── animate.py                    # Animation des barres
│
└── utils/
    ├── __init__.py
    └── profiler.py                   # Mesure du temps, affichage du lagrangien, debug
