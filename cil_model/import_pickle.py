import pickle

# Ouvre le fichier en mode lecture binaire
with open('cil_model/rhs_modifie_num_N2.pkl', 'rb') as fichier:
    donnees = pickle.load(fichier)

# Affiche les données
print(donnees)

