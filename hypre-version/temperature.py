#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Génération de la grille spatiale : N+1 points (N = 100)
x = np.linspace(0, 1, 101)

# Valeurs de u(x) extraites de la sortie du programme
u = np.array([
    1.786432, 1.786432, 1.785621, 1.783980, 1.781472, 1.778040, 1.773607, 1.768074,
    1.761318, 1.753192, 1.743519, 1.732093, 1.718674, 1.702988, 1.684721, 1.663519,
    1.638984, 1.610676, 1.578110, 1.540759, 1.498057, 1.459406, 1.424308, 1.392343,
    1.363156, 1.336442, 1.311939, 1.289423, 1.268695, 1.249585, 1.231940, 1.215627,
    1.200529, 1.186539, 1.173565, 1.161521, 1.150333, 1.139931, 1.130256, 1.121251,
    1.112865, 1.105053, 1.097773, 1.090986, 1.084656, 1.078753, 1.073246, 1.068108,
    1.063314, 1.058841, 1.054667, 1.050772, 1.047138, 1.043748, 1.040586, 1.037637,
    1.034887, 1.032323, 1.029934, 1.027709, 1.025635, 1.023705, 1.021908, 1.020236,
    1.018682, 1.017236, 1.015893, 1.014645, 1.013487, 1.012412, 1.011415, 1.010490,
    1.009634, 1.008840, 1.008105, 1.007425, 1.006796, 1.006214, 1.005676, 1.005179,
    1.004719, 1.004294, 1.003902, 1.003539, 1.003203, 1.002892, 1.002604, 1.002337,
    1.002089, 1.001858, 1.001642, 1.001440, 1.001251, 1.001072, 1.000902, 1.000740,
    1.000584, 1.000434, 1.000287, 1.000143, 1.000000
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

