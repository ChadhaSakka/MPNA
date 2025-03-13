#ifndef DIFFUSION_H
#define DIFFUSION_H

typedef struct {
    double kappa0;
    double sigma;
    double beta;
    double q;  // Exposant dans κ(u) = κ₀ * u^q (déclaré en double pour permettre des valeurs fractionnaires)
} parameters_t;

/**
 * @brief Calcule la valeur de κ(u) = κ₀ * u^q.
 *
 * @param u Valeur de u.
 * @param params Paramètres du problème.
 * @return double Valeur de κ(u).
 */
double compute_kappa(double u, const parameters_t *params);

/**
 * @brief Définit la fonction source Q(x) = β * Heaviside(0.2 - x).
 *
 * @param x Position.
 * @param params Paramètres du problème.
 * @return double Valeur de Q(x).
 */
double source_Q(double x, const parameters_t *params);

/**
 * @brief Initialise la solution u sur la grille.
 *
 * @param u Vecteur de solution (taille N+1).
 * @param N Nombre de points (dx = 1/N).
 * @param params Paramètres du problème.
 * @param dx Pas d'espace.
 */
void discretize_diffusion(double *u, int N, const parameters_t *params, double dx);

/**
 * @brief Met à jour la solution u par une itération du schéma explicite.
 *
 * La mise à jour est effectuée pour les indices internes (i = 1,...,N-1) et
 * les conditions aux limites sont appliquées ensuite : Neumann en i=0 et Dirichlet en i=N.
 *
 * @param u Vecteur de solution courant (taille N+1).
 * @param u_new Vecteur de solution mis à jour.
 * @param N Nombre de points.
 * @param dt Pas de temps.
 * @param params Paramètres du problème.
 * @param dx Pas d'espace.
 */
void update_solution(double *u, double *u_new, int N, double dt, const parameters_t *params, double dx);

#endif // DIFFUSION_H

