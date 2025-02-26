#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_linear_implicit.h"
#include "mesh.h"       // Par exemple, pour computeDx()
#include "utils.h"      // Pour triDiagSolve()

/* 
 * Calcule la conductivité : κ(u) = kappa0 * u^qExponent
 */
static inline double kappa(double kappa0, double qExponent, double u_val) {
    return kappa0 * pow(u_val, qExponent);
}

/* 
 * Calcule la source Q(x) = β * H(δ - x)
 * Ici, δ est fixé à 0.2.
 */
static inline double Q_source(double beta, double x) {
    return (x < 0.2) ? beta : 0.0;
}

void solveDiffusionLinearImplicit(int N, double dt, double finalTime,
                                  double kappa0, double sigma, double beta,
                                  double qExponent,
                                  double *u)
{
    // Domaine : x ∈ [0,1] avec dx = 1/N
    double dx = 1.0 / (double) N;

    // Sauvegarder la solution du pas précédent
    double *u_old = (double*) malloc((N+1) * sizeof(double));
    if(u_old == NULL) { perror("malloc"); exit(EXIT_FAILURE); }

    double time = 0.0;
    while (time < finalTime) {
        if (time + dt > finalTime) {
            dt = finalTime - time;
        }

        // Sauvegarder la solution courante dans u_old
        for (int i = 0; i <= N; i++) {
            u_old[i] = u[i];
        }

        // Allouer les tableaux pour le système tridiagonal de taille N
        double *lower = (double*) malloc(N * sizeof(double));  // indices 1..N-1 utilisés
        double *diag  = (double*) malloc(N * sizeof(double));  // indices 0..N-1
        double *upper = (double*) malloc(N * sizeof(double));  // indices 0..N-2 utilisés
        double *rhs   = (double*) malloc(N * sizeof(double));  // indices 0..N-1
        if (!lower || !diag || !upper || !rhs) { perror("malloc"); exit(EXIT_FAILURE); }

        // Construction du système pour i = 0, …, N–1.
        for (int i = 0; i < N; i++) {
            double x = i * dx;  // position

            // Calcul des valeurs aux interfaces : 
            double u_left, u_right;
            if (i == 0)
                u_left = u_old[1];  // condition de Neumann : u[-1] = u[1]
            else
                u_left = u_old[i-1];
            if (i == N-1)
                u_right = u_old[i+1];  // u[N] est imposé (Dirichlet, u[N]=1)
            else
                u_right = u_old[i+1];

            double kappa_imh = kappa(kappa0, qExponent, 0.5 * (u_old[i] + u_left));   // interface i-1/2
            double kappa_iph = kappa(kappa0, qExponent, 0.5 * (u_old[i] + u_right));  // interface i+1/2

            if (i == 0) {
                // Pour i = 0, condition de Neumann par symétrie
                diag[i]  = 1.0/dt + (2.0/(dx*dx)) * kappa_iph;
                upper[i] = - (2.0/(dx*dx)) * kappa_iph;
                // lower[0] n'est pas utilisé
            } else {
                lower[i] = - (1.0/(dx*dx)) * kappa_imh;
                if (i < N-1) {
                    upper[i] = - (1.0/(dx*dx)) * kappa_iph;
                }
                diag[i] = 1.0/dt + (1.0/(dx*dx)) * (kappa_imh + kappa_iph);
            }

            // Construction du second membre (rhs)
            // Terme d'accumulation implicite, radiation explicitement évaluée à u_old, source Q(x)
            rhs[i] = (1.0/dt) * u_old[i] - sigma * (pow(u_old[i], 4) - 1.0) + Q_source(beta, x);
        }
        // Remarque : l'indice i = N correspond à la condition de Dirichlet u[N]=1.

        // Résoudre le système tridiagonal (taille N)
        triDiagSolve(lower, diag, upper, rhs, N);

        // Mettre à jour la solution pour i = 0, …, N–1, et imposer u[N] = 1.
        for (int i = 0; i < N; i++) {
            u[i] = rhs[i];
        }
        u[N] = 1.0;

        time += dt;

        free(lower);
        free(diag);
        free(upper);
        free(rhs);
    }

    free(u_old);
}

