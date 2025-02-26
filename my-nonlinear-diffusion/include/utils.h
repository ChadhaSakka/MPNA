#ifndef UTILS_H
#define UTILS_H

// A place for shared routines, e.g. a tridiagonal solver
// or boundary condition helpers

// Solve the tridiagonal system A x = rhs
// A is defined by the three diagonals:
//    diagLower[i], diagMain[i], diagUpper[i]
// The solution is stored in x[].
void triDiagSolve(const double *diagLower,
                  const double *diagMain,
                  const double *diagUpper,
                  double *x,  // on input: RHS, on output: solution
                  int n);

#endif

