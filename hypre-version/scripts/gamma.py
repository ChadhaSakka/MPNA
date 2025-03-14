import matplotlib.pyplot as plt
import numpy as np

# Valeurs de gamma
gamma_values = [0.05, 0.10, 0.20]
error_files = [f"erreur_gamma_{gamma:.2f}.txt" for gamma in gamma_values]

# Configuration du style
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(10, 6))

# Couleurs personnalisées
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

# Tracer les courbes
for file, gamma, color in zip(error_files, gamma_values, colors):
    data = np.loadtxt(file, delimiter=',', skiprows=1)
    iterations = data[:, 0]
    errors = data[:, 1]
    plt.semilogy(iterations, errors, label=f'$\gamma = {gamma}$', color=color, lw=2)

# Personnalisation
plt.xlabel("Nombre d'itérations", fontsize=12)
plt.ylabel('Erreur (échelle logarithmique)', fontsize=12)
plt.title('Influence de $\gamma$ sur la convergence (Méthode explicite)', fontsize=14, pad=10)
plt.legend(fontsize=10, frameon=True, edgecolor='black')
plt.grid(True, which="both", linestyle='--', alpha=0.7)
plt.tight_layout()

# Sauvegarde
plt.savefig('gamma_influence.png', dpi=300, bbox_inches='tight')
plt.close()

print("Plot généré : gamma_influence.png")
