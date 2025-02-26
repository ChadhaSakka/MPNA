#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_linear_implicit.h"
#include "diffusion_newton.h"

int main(int argc, char **argv)
{
    // Example usage:
    // 1) read input parameters from CLI or data file
    int N = 50;
    double dt = 0.01;
    double finalTime = 1.0;
    double kappa0 = 0.01;
    double sigma  = 0.1;
    double beta   = 1.0;
    double qExponent = 0.5;
    
    // Allocate solution array (N+1 points)
    double *u = (double*) calloc(N+1, sizeof(double));
    // Initialize with boundary condition: u[N] = 1
    u[N] = 1.0;

    printf("Solving with a linearized implicit scheme...\n");
    solveDiffusionLinearImplicit(N, dt, finalTime, kappa0, sigma, beta, qExponent, u);

    // Print final solution or partial info
    for(int i=0; i<=N; i+=N/5){
        printf("i=%d, u=%.4f\n", i, u[i]);
    }

    // Alternatively, solve with Newton (steady-state)
    // re-init the solution array as an initial guess
    for (int i = 0; i < N+1; i++) {
        u[i] = 1.0; // e.g. initial guess
    }
    int maxNewtonIter = 50;
    double newtonTol  = 1e-8;

    printf("\nSolving with Newton's method for steady state...\n");
    solveDiffusionNewton(N, kappa0, sigma, beta, qExponent,
                         u, maxNewtonIter, newtonTol);

    // Print final solution
    for(int i=0; i<=N; i+=N/5){
        printf("i=%d, u=%.4f\n", i, u[i]);
    }

    free(u);
    return 0;
}

