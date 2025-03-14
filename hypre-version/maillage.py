import matplotlib.pyplot as plt
import numpy as np

# Valeurs de N
N_values = [50, 100, 200]
profile_files = [f"profil_N_{N}.txt" for N in N_values]

# Configuration du style
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(10, 6))

# Couleurs personnalisées
colors = ['#9467bd', '#d62728', '#8c564b']

# Tracer les courbes
for file, N, color in zip(profile_files, N_values, colors):
    data = np.loadtxt(file, delimiter=',', skiprows=1)
    x = data[:, 0]
    u = data[:, 1]
    plt.plot(x, u, label=f'$N = {N}$', color=color, lw=2, marker='o', markersize=4, markevery=N//10)

# Personnalisation
plt.xlabel('Position $x$', fontsize=12)
plt.ylabel('Température $u(x)$', fontsize=12)
plt.title('Influence du maillage sur le profil (Méthode de Newton)', fontsize=14, pad=10)
plt.legend(fontsize=10, frameon=True, edgecolor='black')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()

# Sauvegarde
plt.savefig('mesh_influence.png', dpi=300, bbox_inches='tight')
plt.close()

print("Plot généré : mesh_influence.png")
