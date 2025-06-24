# models/base_system.py

from abc import ABC, abstractmethod

class BaseSystem(ABC):
    """
    Classe de base abstraite pour un système de N barres rigides.
    Ne contient que le strict nécessaire commun.
    """

    def __init__(self, N):
        """
        Initialise le nombre de barres.

        :param N: Nombre de barres
        """
        self.N = N

    @abstractmethod
    def solve(self, t_span, t_eval=None, method='RK45'):
        """
        Méthode abstraite : doit être implémentée dans chaque sous-classe.
        Permet de résoudre la dynamique du système.

        :param t_span: Intervalle de temps (t0, tf)
        :param t_eval: Points de temps pour l'échantillonnage
        :param method: Méthode d'intégration (ex: 'RK45', 'DOP853', etc.)
        """
        pass
