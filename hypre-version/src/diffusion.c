#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion.h"

/**
 * @brief Calcule le coefficient de diffusion κ(u) = κ₀ * u^q.
 */
double compute_kappa(double u, const parameters_t *params) {
    return params->kappa0 * pow(u, params->q);
}

/**
 * @brief Calcule le terme source Q(x) = β * Heaviside(0.2 - x).
 */
double source_Q(double x, const parameters_t *params) {
    return (x < 0.2) ? params->beta : 0.0;
}

/**
 * @brief Initialise la solution sur la grille : u = 1 partout.
 */
void discretize_diffusion(double *u, int N, const parameters_t *params, double dx) {
    for (int i = 0; i <= N; i++) {
        u[i] = 1.0;
    }
    // Les conditions aux limites seront appliquées lors des itérations :
    // Neumann en i=0 et Dirichlet en i=N.
}

/**
 * @brief Met à jour la solution u en utilisant le schéma explicite.
 *
 * La mise à jour est effectuée pour i = 1,...,N-1. Les conditions aux limites sont
 * ensuite imposées explicitement :
 * - Neumann en i = 0 : u_new[0] = u_new[1]
 * - Dirichlet en i = N : u_new[N] = 1.0
 */
void update_solution(double *u, double *u_new, int N, double dt, const parameters_t *params, double dx) {
    for (int i = 1; i < N; i++) {
        double u_ip = u[i+1];
        double u_im = u[i-1];  // Pour i>=1, u[i-1] est défini
        // Calcul des coefficients de diffusion (moyenne sur les points voisins)
        double kappa_ip = compute_kappa((u[i+1] + u[i]) / 2.0, params);
        double kappa_im = compute_kappa((u[i-1] + u[i]) / 2.0, params);
        // Discrétisation du terme de diffusion
        double diffusion_term = (kappa_ip * (u_ip - u[i]) - kappa_im * (u[i] - u_im)) / (dx * dx);
        // Terme de rayonnement
        double radiation_term = params->sigma * (pow(u[i], 4) - 1.0);
        // Terme source
        double x = i * dx;
        double Q = source_Q(x, params);
        // Mise à jour de la solution
        u_new[i] = u[i] + dt * (diffusion_term - radiation_term + Q);
    }
    // Application explicite des conditions aux limites
    u_new[0] = u_new[1];  // Neumann en i=0
    u_new[N] = 1.0;       // Dirichlet en i=N
}

