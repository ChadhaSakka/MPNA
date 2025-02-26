#ifndef DIFFUSION_LINEAR_IMPLICIT_H
#define DIFFUSION_LINEAR_IMPLICIT_H

/* 
 * Solve the nonlinear diffusion problem using
 * a "linearized implicit" time-integration approach.
 *
 *  - N: number of space discretization points
 *  - dt: time step
 *  - finalTime: final simulation time
 *  - kappa0, sigma, beta: physical parameters
 *  - qExponent: exponent in kappa(u) = kappa0 * u^qExponent
 *  - u[]: array to store the solution
 */
void solveDiffusionLinearImplicit(int N, double dt, double finalTime,
                                  double kappa0, double sigma, double beta,
                                  double qExponent,
                                  double *u);

#endif

