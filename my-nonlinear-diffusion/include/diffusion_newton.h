#ifndef DIFFUSION_NEWTON_H
#define DIFFUSION_NEWTON_H

/* 
 * Solve the nonlinear diffusion problem using
 * Newton’s method on the steady-state equations.
 *
 *  - N: number of space discretization points
 *  - kappa0, sigma, beta, qExponent: physical parameters
 *  - u[]: array to store the solution (initial guess in input, final solution in output)
 *  - maxNewtonIter: maximum number of Newton iterations
 *  - newtonTol: tolerance for stopping
 */
void solveDiffusionNewton(int N,
                          double kappa0, double sigma, double beta,
                          double qExponent,
                          double *u,
                          int maxNewtonIter, double newtonTol);

#endif

