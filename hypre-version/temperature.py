#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Génération de la grille spatiale : N+1 points (N = 100)
x = np.linspace(0, 1, 101)

# Valeurs de u(x) extraites de la nouvelle sortie du programme
u = np.array([
    1.771811, 1.771811, 1.770964, 1.769256, 1.766656, 1.763120, 1.758587, 1.752979,
    1.746198, 1.738130, 1.728635, 1.717552, 1.704693, 1.689840, 1.672743, 1.653116,
    1.630632, 1.604919, 1.575551, 1.542045, 1.503850, 1.468599, 1.436009, 1.405833,
    1.377855, 1.351882, 1.327744, 1.305288, 1.284379, 1.264896, 1.246729, 1.229778,
    1.213954, 1.199177, 1.185372, 1.172471, 1.160413, 1.149142, 1.138605, 1.128755,
    1.119546, 1.110940, 1.102897, 1.095382, 1.088363, 1.081809, 1.075692, 1.069985,
    1.064663, 1.059702, 1.055081, 1.050779, 1.046777, 1.043055, 1.039597, 1.036387,
    1.033408, 1.030647, 1.028089, 1.025722, 1.023533, 1.021511, 1.019644, 1.017922,
    1.016336, 1.014876, 1.013533, 1.012299, 1.011167, 1.010129, 1.009178, 1.008307,
    1.007512, 1.006785, 1.006122, 1.005518, 1.004968, 1.004467, 1.004012, 1.003599,
    1.003224, 1.002884, 1.002576, 1.002297, 1.002044, 1.001815, 1.001608, 1.001420,
    1.001250, 1.001096, 1.000955, 1.000827, 1.000709, 1.000601, 1.000501, 1.000407,
    1.000319, 1.000236, 1.000155, 1.000077, 1.000000
])

# Style et amélioration esthétique
plt.style.use('ggplot')
plt.figure(figsize=(10, 6))
plt.plot(x, u, marker='o', markersize=5, linewidth=2, label='Température')
plt.xlabel("x", fontsize=14)
plt.ylabel("u(x)", fontsize=14)
plt.title("Profil de température $u(x)$", fontsize=16, fontweight='bold')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(fontsize=12)

# Annotation sur le point caractéristique x = 0.2
idx = np.argmin(np.abs(x - 0.2))
plt.annotate(f"u(0.2) = {u[idx]:.3f}",
             xy=(x[idx], u[idx]), xytext=(0.3, u[idx] + 0.1),
             arrowprops=dict(facecolor='black', shrink=0.05),
             fontsize=12)

plt.tight_layout()
plt.savefig('resultats_improved.png', dpi=300)
plt.show()

