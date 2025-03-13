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
    params.q = 0.5; // sqrt(u)

    int N = 10;
    double dx = 1.0 / N;
    double *u = alloc_vector(N+1);

    discretize_diffusion(u, N, &params, dx);
    printf("Test de discretize_diffusion:\n");
    print_vector(u, N+1);

    free_vector(u);
    return 0;
}

