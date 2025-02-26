#include <stdio.h>
#include <stdlib.h>
#include "diffusion_linear_implicit.h"

int main(void) {
    int N = 10;
    double dt = 0.01;
    double finalTime = 0.1;
    // Paramètres physiques (exemple)
    double kappa0 = 0.01, sigma = 0.1, beta = 1.0, qExponent = 0.5;

    double *u = (double*) calloc(N+1, sizeof(double));
    if(u == NULL) { perror("calloc"); exit(EXIT_FAILURE); }
    // Initialisation : on peut commencer par u=0 pour i=0,...,N-1 et u[N]=1
    for (int i = 0; i < N; i++) {
        u[i] = 0.0;
    }
    u[N] = 1.0;

    solveDiffusionLinearImplicit(N, dt, finalTime, kappa0, sigma, beta, qExponent, u);

    for (int i = 0; i <= N; i++) {
        printf("i=%d, u=%.6f\n", i, u[i]);
    }
    
    free(u);
    return 0;
}

