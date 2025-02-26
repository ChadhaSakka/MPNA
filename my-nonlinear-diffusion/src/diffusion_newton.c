#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_newton.h"
#include "mesh.h"
#include "utils.h"

void solveDiffusionNewton(int N,
                          double kappa0, double sigma, double beta,
                          double qExponent,
                          double *u,
                          int maxNewtonIter, double newtonTol)
{
    double dx = 1.0 / (double)N;

    // For each Newton iteration:
    for (int iter = 0; iter < maxNewtonIter; iter++) {
        
        // 1) Build the Jacobian (tridiagonal) and the residual vector:
        //    J * deltaU = -F(u)
        //    where F(u) = 0 are the discrete equations.

        // 2) Solve the linear system for deltaU:
        //    triDiagSolve(a, b, c, rhs, deltaU)

        // 3) Update: u <- u + deltaU

        // 4) Check convergence: if ||deltaU|| < newtonTol, break
    }
}

