#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "diffusion.h"
#include "solveur.h"
#include "utils.h"

int main(int argc, char **argv) {	
    MPI_Init(&argc, &argv);
    
    // Choix des paramètres (Exemple 1)
    parameters_t params;
    params.kappa0 = 0.01;
    params.sigma = 0.1;
    params.beta = 1.0;
    params.q = 0.5;  // Pour κ(u) = κ₀ * u^(q) (sqrt(u) si q=0.5)

    int N = 100;         // Nombre de points (modifiable entre 50 et 10000)
    double dx = 1.0 / N; // Pas spatial
    double gamma = 0.1;  // Facteur pour dt (peut être 0.1, 1, ou 10)

    // Allocation des vecteurs solution
    double *u = alloc_vector(N+1);
    double *u_new = alloc_vector(N+1);

    // Initialisation de la solution
    discretize_diffusion(u, N, &params, dx);

    // Affichage des paramètres et informations initiales
    printf("=== Résolution de l'équation de diffusion non linéaire ===\n");
    printf("Paramètres:\n");
    printf("  κ₀ = %g\n", params.kappa0);
    printf("  σ  = %g\n", params.sigma);
    printf("  β  = %g\n", params.beta);
    printf("  q  = %g\n", params.q);
    printf("Grille: N = %d, dx = %g\n", N, dx);
    printf("Facteur dt (gamma) = %g\n\n", gamma);

    // Entête de la table d'itération
    printf("%-10s %-15s %-15s\n", "Itération", "dt", "Différence");

    // Boucle d'itération jusqu'à convergence (stationnaire)
    double tol = 1e-6;
    int max_iter = 10000;
    int iter = 0;
    double diff = 1.0;
    while(diff > tol && iter < max_iter) {
        // Calcul du pas en temps via la fonction dédiée
        double dt = compute_dt(u, N, &params, dx, gamma);
        
        // Mise à jour de la solution
        update_solution(u, u_new, N, dt, &params, dx);
        
        // Calcul de la différence pour vérifier la convergence et mise à jour de u
        diff = 0.0;
        for (int i = 0; i <= N; i++) {
            diff += fabs(u_new[i] - u[i]);
            u[i] = u_new[i];
        }
        iter++;
        // Affichage toutes les 100 itérations
        if(iter % 100 == 0)
            printf("%-10d %-15.6g %-15.6e\n", iter, dt, diff);
    }
    printf("\nConvergence atteinte en %d itérations, diff = %e\n\n", iter, diff);
    
    // Affichage de la solution finale sous forme de table
    printf("Solution finale u(x):\n");
    printf("%-10s %-15s\n", "x", "u(x)");
    for (int i = 0; i <= N; i++) {
        double x_value = i * dx;
        printf("%-10.4f %-15.6f\n", x_value, u[i]);
    }
    
    // Libération de la mémoire
    free_vector(u);
    free_vector(u_new);
    
    MPI_Finalize();
    return 0;
}

