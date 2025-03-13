#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Données extraites de la sortie du programme (toutes les 100 itérations)
iterations = np.arange(100, 10001, 100)
diff_values = np.array([
    9.778997e-03, 9.543144e-03, 9.286071e-03, 9.009255e-03, 8.714685e-03, 8.404973e-03,
    8.083109e-03, 7.752237e-03, 7.415502e-03, 7.075949e-03, 6.736452e-03, 6.399651e-03,
    6.067917e-03, 5.743319e-03, 5.427611e-03, 5.122233e-03, 4.828316e-03, 4.546703e-03,
    4.277970e-03, 4.022457e-03, 3.780293e-03, 3.551436e-03, 3.335691e-03, 3.132749e-03,
    2.942207e-03, 2.763590e-03, 2.596376e-03, 2.440011e-03, 2.293919e-03, 2.157522e-03,
    2.030242e-03, 1.911512e-03, 1.800783e-03, 1.697523e-03, 1.601227e-03, 1.511413e-03,
    1.427626e-03, 1.349436e-03, 1.276442e-03, 1.208269e-03, 1.144564e-03, 1.085004e-03,
    1.029283e-03, 9.771229e-04, 9.282618e-04, 8.824596e-04, 8.394939e-04, 7.991593e-04,
    7.612661e-04, 7.256395e-04, 6.921179e-04, 6.605524e-04, 6.308057e-04, 6.027512e-04,
    5.762720e-04, 5.512602e-04, 5.276163e-04, 5.052485e-04, 4.840719e-04, 4.640082e-04,
    4.449849e-04, 4.269351e-04, 4.097967e-04, 3.935123e-04, 3.780289e-04, 3.632973e-04,
    3.492716e-04, 3.359096e-04, 3.231719e-04, 3.110220e-04, 2.994258e-04, 2.883517e-04,
    2.777703e-04, 2.676542e-04, 2.579778e-04, 2.487172e-04, 2.398501e-04, 2.313557e-04,
    2.232145e-04, 2.154083e-04, 2.079201e-04, 2.007338e-04, 1.938345e-04, 1.872080e-04,
    1.808412e-04, 1.747216e-04, 1.688376e-04, 1.631781e-04, 1.577327e-04, 1.524917e-04,
    1.474458e-04, 1.425864e-04, 1.379053e-04, 1.333946e-04, 1.290470e-04, 1.248556e-04,
    1.208137e-04, 1.169151e-04, 1.131538e-04, 1.095243e-04
])

# Amélioration esthétique avec un style moderne
plt.style.use('ggplot')
plt.figure(figsize=(10, 6))
plt.semilogy(iterations, diff_values, marker='o', markersize=5, linewidth=2, label='Différence')
# Ligne horizontale pour la tolérance finale
plt.axhline(y=1e-4, color='black', linestyle='--', linewidth=1.5, label='Tolérance (1e-4)')

plt.xlabel("Nombre d'itérations", fontsize=14)
plt.ylabel("Différence (échelle logarithmique)", fontsize=14)
plt.title("Convergence de la solution", fontsize=16, fontweight='bold')
plt.legend(fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('convergence_improved.png', dpi=300)
plt.show()

