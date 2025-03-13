#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Données extraites de la nouvelle sortie du programme (toutes les 100 itérations)
iterations = np.arange(100, 10001, 100)
diff_values = np.array([
    9.553546e-03, 9.130764e-03, 8.724190e-03, 8.332150e-03, 7.953708e-03, 7.588520e-03,
    7.236493e-03, 6.897580e-03, 6.571712e-03, 6.258775e-03, 5.958618e-03, 5.671058e-03,
    5.395879e-03, 5.132840e-03, 4.881678e-03, 4.642102e-03, 4.413805e-03, 4.196456e-03,
    3.989713e-03, 3.793215e-03, 3.606593e-03, 3.429471e-03, 3.261466e-03, 3.102193e-03,
    2.951271e-03, 2.808317e-03, 2.672957e-03, 2.544823e-03, 2.423556e-03, 2.308805e-03,
    2.200232e-03, 2.097511e-03, 2.000327e-03, 1.908379e-03, 1.821377e-03, 1.739048e-03,
    1.661127e-03, 1.587366e-03, 1.517527e-03, 1.451388e-03, 1.388734e-03, 1.329365e-03,
    1.273091e-03, 1.219734e-03, 1.169125e-03, 1.121104e-03, 1.075523e-03, 1.032240e-03,
    9.911237e-04, 9.520486e-04, 9.148982e-04, 8.795625e-04, 8.459383e-04, 8.139287e-04,
    7.834425e-04, 7.543944e-04, 7.267040e-04, 7.002959e-04, 6.750993e-04, 6.510476e-04,
    6.280786e-04, 6.061334e-04, 5.851569e-04, 5.650974e-04, 5.459062e-04, 5.275374e-04,
    5.099481e-04, 4.930977e-04, 4.769482e-04, 4.614638e-04, 4.466107e-04, 4.323571e-04,
    4.186731e-04, 4.055306e-04, 3.929029e-04, 3.807650e-04, 3.690933e-04, 3.578655e-04,
    3.470605e-04, 3.366584e-04, 3.266406e-04, 3.169892e-04, 3.076876e-04, 2.987198e-04,
    2.900709e-04, 2.817266e-04, 2.736736e-04, 2.658992e-04, 2.583912e-04, 2.511383e-04,
    2.441297e-04, 2.373550e-04, 2.308045e-04, 2.244689e-04, 2.183395e-04, 2.124079e-04,
    2.066662e-04, 2.011069e-04, 1.957227e-04, 1.905067e-04
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

