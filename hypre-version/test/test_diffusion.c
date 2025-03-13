#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion.h"
#include "utils.h"

int main() {
    parameters_t params;
    params.kappa0 = 0.01;
    params.sigma = 0.1;
    params.beta = 1.0;
    params.q = 0.5; // Pour κ(u) = κ₀ * u^(q)

    int N = 10;
    double dx = 1.0 / N;
    double *u = alloc_vector(N+1);

    // Initialisation de la solution (u = 1 partout)
    discretize_diffusion(u, N, &params, dx);  // Correction : &params au lieu de ¶ms
    printf("Test de discretize_diffusion:\n");
    print_vector(u, N+1);

    free_vector(u);
    return 0;
}
