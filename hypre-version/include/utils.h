#ifndef UTILS_H
#define UTILS_H

#include "diffusion.h"  // Ajouté pour la définition de parameters_t et compute_kappa

/**
 * @brief Alloue dynamiquement un vecteur de taille n et l'initialise à zéro.
 *
 * @param n Taille du vecteur.
 * @return double* Pointeur sur le vecteur alloué.
 */
double* alloc_vector(int n);

/**
 * @brief Libère la mémoire allouée pour un vecteur.
 *
 * @param v Pointeur sur le vecteur.
 */
void free_vector(double *v);

/**
 * @brief Affiche un vecteur de taille n.
 *
 * @param v Vecteur à afficher.
 * @param n Taille du vecteur.
 */
void print_vector(const double *v, int n);

/**
 * @brief Calcule le pas de temps dt en fonction de la solution u.
 *
 * dt est calculé à partir de la condition de stabilité :
 * dt < 2 / (4σ u_max^3 + 4κ(u_max)/(dx^2))
 *
 * @param u      Vecteur solution (taille N+1).
 * @param N      Nombre de points.
 * @param params Paramètres du problème.
 * @param dx     Pas spatial.
 * @param gamma  Facteur multiplicatif.
 * @return double dt calculé.
 */
double compute_dt(const double *u, int N, const parameters_t *params, double dx, double gamma);

#endif // UTILS_H

