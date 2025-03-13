#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion.h"  // Inclusion pour la définition de parameters_t et compute_kappa
#include "utils.h"

/**
 * @brief Alloue et initialise un vecteur de taille n.
 *
 * @param n Taille du vecteur.
 * @return double* Pointeur vers le vecteur alloué.
 */
double* alloc_vector(int n) {
    double *v = (double*) malloc(n * sizeof(double));
    if (!v) {
        fprintf(stderr, "Erreur d'allocation pour un vecteur de taille %d\n", n);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++) {
        v[i] = 0.0;
    }
    return v;
}

/**
 * @brief Libère la mémoire allouée pour un vecteur.
 */
void free_vector(double *v) {
    free(v);
}

/**
 * @brief Affiche le contenu d'un vecteur.
 *
 * @param v Pointeur sur le vecteur.
 * @param n Taille du vecteur.
 */
void print_vector(const double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%f ", v[i]);
    }
    printf("\n");
}

/**
 * @brief Calcule le pas de temps dt en fonction de la solution u.
 *
 * dt est calculé à partir de la condition de stabilité :
 * dt < 2 / (4σ u_max^3 + 4κ(u_max)/(dx^2))
 *
 * @param u      Vecteur solution (taille N+1).
 * @param N      Nombre de points.
 * @param params Paramètres physiques.
 * @param dx     Pas spatial.
 * @param gamma  Facteur multiplicatif.
 * @return double dt calculé.
 */
double compute_dt(const double *u, int N, const parameters_t *params, double dx, double gamma) {
    double u_max = 0.0;
    for (int i = 0; i <= N; i++) {
        if (fabs(u[i]) > u_max) {
            u_max = fabs(u[i]);
        }
    }
    return gamma * 2.0 / (4.0 * params->sigma * pow(u_max, 3) +
                            4.0 * compute_kappa(u_max, params) / (dx * dx));
}

