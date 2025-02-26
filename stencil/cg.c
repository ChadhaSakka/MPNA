#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*-------------------------------------------------------------
 * Helper functions for vector/matrix operations
 *-------------------------------------------------------------*/

/* Matrix-vector product: y = A*x
 * A is in row-major order, dimension n x n
 * x, y are vectors of length n
 */
void matVecProduct(const double *A, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i*n + j] * x[j];
        }
        y[i] = sum;
    }
}

/* Dot product of two vectors x, y of length n */
double dotProduct(const double *x, const double *y, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}

/* Scaled vector addition: y = y + alpha * x
 * (a.k.a. AXPY in BLAS)
 */
void axpy(double alpha, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

/* Euclidian norm of a vector */
double vecNorm(const double *x, int n) {
    return sqrt(dotProduct(x, x, n));
}

/*-------------------------------------------------------------
 * Conjugate Gradient solver
 *-------------------------------------------------------------*/

/* Solves A*x = b using Conjugate Gradient
 *  - A: pointer to matrix (size n*n)
 *  - b: pointer to RHS vector (size n)
 *  - x: pointer to initial guess + solution (size n)
 *  - n: dimension
 *  - tol: tolerance for stopping condition
 *  - maxIter: maximum number of iterations
 * 
 * Returns the number of iterations used, or -1 if it fails.
 */
int conjugateGradient(const double *A, const double *b,
                      double *x, int n,
                      double tol, int maxIter)
{
    double *r = (double*) malloc(n * sizeof(double));
    double *p = (double*) malloc(n * sizeof(double));
    double *Ap = (double*) malloc(n * sizeof(double));

    // r0 = b - A*x0
    matVecProduct(A, x, r, n);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    double rDotr = dotProduct(r, r, n);
    double rNorm = sqrt(rDotr);
    double rNorm0 = rNorm;  // for relative convergence check

    int k = 0;
    for (; k < maxIter; k++) {
        // Check convergence: relative residual
        if (rNorm / (rNorm0 + 1e-14) < tol) {
            break;
        }

        // Ap = A * p
        matVecProduct(A, p, Ap, n);

        // alpha = (r^T r) / (p^T A p)
        double alpha = rDotr / dotProduct(p, Ap, n);

        // x_{k+1} = x_k + alpha * p
        axpy(alpha, p, x, n);

        // r_{k+1} = r_k - alpha * Ap
        axpy(-alpha, Ap, r, n);

        double rDotr_new = dotProduct(r, r, n);

        // beta = (r_{k+1}^T r_{k+1}) / (r_k^T r_k)
        double beta = rDotr_new / rDotr;
        rDotr = rDotr_new;

        // p_{k+1} = r_{k+1} + beta * p_k
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        rNorm = sqrt(rDotr);
    }

    free(r);
    free(p);
    free(Ap);

    // If we exited by reaching maxIter, you can check if we converged
    if (k == maxIter) {
        // might not be a failure necessarily, but let's say so
        return -1;
    }

    return k;
}

/*-------------------------------------------------------------
 * Example main (demo usage)
 *-------------------------------------------------------------*/
int main(void) {
    int n = 4; // dimension
    double A[16] = {
        4, 1, 0, 0,
        1, 3, 1, 0,
        0, 1, 3, 1,
        0, 0, 1, 2
    };
    double b[4] = { 1, 2, 2, 1 };

    double x[4] = { 0.0, 0.0, 0.0, 0.0 }; // initial guess

    double tol = 1e-8;
    int maxIter = 1000;

    int iters = conjugateGradient(A, b, x, n, tol, maxIter);
    if (iters < 0) {
        printf("CG did not converge within %d iterations.\n", maxIter);
    } else {
        printf("CG converged in %d iterations.\n", iters);
        printf("Solution x = [");
        for (int i = 0; i < n; i++) {
            printf(" %g", x[i]);
        }
        printf(" ]\n");
    }

    return 0;
}

