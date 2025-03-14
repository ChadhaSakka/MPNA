import matplotlib.pyplot as plt
import numpy as np

# Charger les données
data = np.loadtxt("solutions.csv", delimiter=',', skiprows=1)
x = data[:, 0]
u_explicit = data[:, 1]
u_implicit = data[:, 2]
u_newton = data[:, 3]

# Données de performance (tirées de vos résultats)
methods = ['Explicite', 'Implicite', 'Newton']
iterations = [10000, 1000, 7]
errors = [1.886821e-06, 6.308360e-05, 1.979606e-07]
times = [0.113882, 0.0756467, 0.0043864]

# Configuration du style
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(10, 6))

# Couleurs et styles
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Bleu, Orange, Vert
markers = ['o', 's', '^']  # Cercle, Carré, Triangle
linestyles = ['-', '--', '-.']  # Solide, Tirets, Pointillés

# Tracer les courbes
for i, (u, method, color, marker, ls) in enumerate(zip([u_explicit, u_implicit, u_newton], methods, colors, markers, linestyles)):
    plt.plot(x, u, label=f'{method} ({iterations[i]} itér., {errors[i]:.1e}, {times[i]:.3f}s)', 
             color=color, lw=2.5, linestyle=ls, marker=marker, markersize=6, markevery=10, alpha=0.85)

# Personnalisation
plt.xlabel('Position $x$', fontsize=14)
plt.ylabel('Température $u(x)$', fontsize=14)
plt.title('Comparaison des méthodes pour $N = 100$, $\gamma = 0.1$', fontsize=16, pad=15)
plt.legend(fontsize=11, loc='upper right', frameon=True, edgecolor='black', bbox_to_anchor=(1.0, 1.0))
plt.grid(True, linestyle='--', alpha=0.5)
plt.xlim(0, 1)
plt.ylim(0.95, 1.85)  # Ajusté pour mieux voir les différences

# Ajustements finaux
plt.tight_layout()

# Sauvegarde
plt.savefig('method_comparison_clear.png', dpi=300, bbox_inches='tight')
plt.close()

print("Plot généré : method_comparison_clear.png")
