#include <stdlib.h>
#include "utils.h"

/*
 * Algorithme de Thomas pour résoudre un système tridiagonal de taille n.
 * - lower : tableau de taille n (indice 0 non utilisé, valides pour i=1..n-1)
 * - diag  : tableau de taille n (indices 0..n-1)
 * - upper : tableau de taille n (indice n-1 non utilisé, valides pour i=0..n-2)
 * - b     : second membre (taille n); la solution sera écrite dans b.
 */
void triDiagSolve(const double *lower, const double *diag, const double *upper, double *b, int n) {
    double *c = (double*) malloc(n * sizeof(double));
    double *d = (double*) malloc(n * sizeof(double));

    d[0] = diag[0];
    if (d[0] == 0.0) { free(c); free(d); return; }
    c[0] = upper[0] / d[0];

    for (int i = 1; i < n; i++) {
        d[i] = diag[i] - lower[i] * c[i-1];
        if (d[i] == 0.0) { free(c); free(d); return; }
        if (i < n - 1) {
            c[i] = upper[i] / d[i];
        }
    }

    // Substitution avant (modification du second membre)
    for (int i = 1; i < n; i++) {
        b[i] = b[i] - lower[i] * b[i-1];
    }
    // Substitution arrière
    b[n-1] = b[n-1] / d[n-1];
    for (int i = n - 2; i >= 0; i--) {
        b[i] = (b[i] - upper[i] * b[i+1]) / d[i];
    }
    
    free(c);
    free(d);
}

