#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion.h"

/**
 * @brief Calcule le coefficient de diffusion κ(u) = κ₀ * u^(q).
 *
 * @param u La valeur de u.
 * @param params Les paramètres physiques (incluant κ₀ et l’exposant q).
 * @return double La valeur de κ(u).
 */
double compute_kappa(double u, const parameters_t *params) {
    return params->kappa0 * pow(u, params->q);
}

/**
 * @brief Calcule le terme source Q(x) = β * Heaviside(0.2 - x).
 *
 * @param x La position spatiale.
 * @param params Les paramètres physiques (incluant β).
 * @return double La valeur de Q(x).
 */
double source_Q(double x, const parameters_t *params) {
    return (x < 0.2) ? params->beta : 0.0;
}

/**
 * @brief Initialise la solution sur la grille : u = 1 partout.
 *
 * @param u Le vecteur de solution (taille N+1).
 * @param N Le nombre de segments (la grille contient N+1 points).
 * @param params Les paramètres du problème (non utilisés ici, mais conservés pour cohérence).
 * @param dx Le pas spatial.
 */
void discretize_diffusion(double *u, int N, const parameters_t *params, double dx) {
    for (int i = 0; i <= N; i++) {
        u[i] = 1.0;
    }
    // Les conditions aux limites seront imposées lors de la résolution.
}

/**
 * @brief Met à jour la solution u en utilisant le schéma explicite.
 *
 * La mise à jour est effectuée pour les noeuds internes (i = 1, …, N-1). 
 * Les conditions aux limites sont ensuite imposées explicitement :
 * - Neumann en i=0 : u_new[0] = u_new[1]
 * - Dirichlet en i=N : u_new[N] = 1.0
 *
 * @param u Le vecteur solution courant.
 * @param u_new Le vecteur solution mis à jour.
 * @param N Le nombre de segments (la grille contient N+1 points).
 * @param dt Le pas de temps.
 * @param params Les paramètres du problème.
 * @param dx Le pas spatial.
 */
void update_solution(double *u, double *u_new, int N, double dt, const parameters_t *params, double dx) {
    for (int i = 1; i < N; i++) {
        double u_ip = u[i+1];
        double u_im = u[i-1];
        // Calcul correct de κ_{i+1/2} et κ_{i-1/2}
        double kappa_ip = 0.5 * (compute_kappa(u[i], params) + compute_kappa(u[i+1], params));
        double kappa_im = 0.5 * (compute_kappa(u[i-1], params) + compute_kappa(u[i], params));
        // Discrétisation du terme de diffusion
        double diffusion_term = (kappa_ip * (u_ip - u[i]) - kappa_im * (u[i] - u_im)) / (dx * dx);
        // Terme de rayonnement
        double radiation_term = params->sigma * (pow(u[i], 4) - 1.0);
        // Terme source
        double x = i * dx;
        double Q = source_Q(x, params);
        // Mise à jour explicite de la solution
        u_new[i] = u[i] + dt * (diffusion_term - radiation_term + Q);
    }
    // Application explicite des conditions aux limites
    u_new[0] = u_new[1];  // Neumann en x=0
    u_new[N] = 1.0;       // Dirichlet en x=1
}
