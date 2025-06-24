import pickle
import os

# Affiche le répertoire courant
print("Répertoire de travail :", os.getcwd())

# Chemin vers le fichier (à adapter)
chemin_fichier = "rhs_modifie_num_N3.pkl"  # Si dans le même dossier

# Lecture du fichier
try:
    with open(chemin_fichier, 'rb') as fichier:
        donnees = pickle.load(fichier)
    print("Données chargées :", donnees)
except FileNotFoundError:
    print("Erreur : Fichier non trouvé. Vérifie le chemin :", chemin_fichier)