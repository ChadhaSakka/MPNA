#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def plot_temperature_profile():
    """
    Trace le profil de température u(x) en utilisant les données extraites de la simulation.
    """
    # Génération de la grille spatiale (N=100 -> 101 points)
    x = np.linspace(0, 1, 101)
    # Valeurs de u(x) extraites de la dernière simulation (à adapter si besoin)
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

def plot_convergence():
    """
    Trace l'évolution de la convergence (la différence entre itérations) en échelle semi-logarithmique.
    """
    # Itérations (toutes les 100 itérations de 100 à 10000)
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

def plot_gamma_comparison():
    """
    Compare la convergence pour différentes valeurs de gamma.
    Ici, nous générons des données synthétiques pour illustrer l'influence de gamma sur le taux de décroissance.
    """
    plt.style.use('ggplot')
    plt.figure(figsize=(10, 6))
    
    iterations = np.linspace(0, 10000, 100)
    gammas = [0.1, 1, 10]
    for gamma in gammas:
        # Données synthétiques : décroissance exponentielle avec un taux influencé par gamma
        error = 0.01 * np.exp(-gamma * iterations / 10000)
        plt.semilogy(iterations, error, label=f'gamma = {gamma}', linewidth=2)
    
    plt.xlabel("Nombre d'itérations", fontsize=14)
    plt.ylabel("Erreur", fontsize=14)
    plt.title("Comparaison de convergence pour différentes valeurs de gamma", fontsize=16, fontweight='bold')
    plt.legend(fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig('gamma_comparison.png', dpi=300)
    plt.show()

def plot_mesh_comparison():
    """
    Compare le profil de température pour différentes tailles de maillage.
    Ici, nous générons des données synthétiques pour illustrer l'influence de N sur la solution.
    """
    plt.style.use('ggplot')
    plt.figure(figsize=(10, 6))
    
    # Différents nombres de points (maillages)
    meshes = [50, 100, 500]
    for N in meshes:
        x = np.linspace(0, 1, N+1)
        # Pour l'exemple, simulons u(x) comme une courbe décroissante avec un léger ajustement lié au maillage
        # La solution tend vers 1 pour x=1 et prend une valeur plus élevée au début.
        u = 1 + 0.77 * np.cos(np.pi * x) * (1 + 1/N)
        plt.plot(x, u, marker='o', markersize=3, linewidth=2, label=f'N = {N}')
    
    plt.xlabel("x", fontsize=14)
    plt.ylabel("u(x)", fontsize=14)
    plt.title("Comparaison des profils de température pour différents maillages", fontsize=16, fontweight='bold')
    plt.legend(fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig('mesh_comparison.png', dpi=300)
    plt.show()

def generate_all_plots():
    plot_temperature_profile()
    plot_convergence()
    plot_gamma_comparison()
    plot_mesh_comparison()

if __name__ == "__main__":
    generate_all_plots()

