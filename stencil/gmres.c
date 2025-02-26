#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*-------------------------------------------------------------
 * Fonctions utilitaires
 *-------------------------------------------------------------*/

// Produit matrice-vecteur : y = A * x
// A est stockée en ordre ligne (row-major) et est de dimension n x n
void matVecProduct(const double *A, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        double somme = 0.0;
        for (int j = 0; j < n; j++) {
            somme += A[i*n + j] * x[j];
        }
        y[i] = somme;
    }
}

// Produit scalaire de deux vecteurs de longueur n
double dotProduct(const double *x, const double *y, int n) {
    double somme = 0.0;
    for (int i = 0; i < n; i++) {
        somme += x[i] * y[i];
    }
    return somme;
}

// Opération y = y + alpha * x (axpy)
void axpy(double alpha, const double *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

// Calcule la norme euclidienne d'un vecteur de longueur n
double vecNorm(const double *x, int n) {
    return sqrt(dotProduct(x, x, n));
}

/*-------------------------------------------------------------
 * Solveur GMRES (Full GMRES)
 *-------------------------------------------------------------*/

/*
 * Résout A*x = b en utilisant le GMRES complet.
 *
 *  - A : matrice de dimension n x n, stockée en ordre ligne.
 *  - b : vecteur second membre de longueur n.
 *  - x : vecteur initial (et solution en sortie) de longueur n.
 *  - n : dimension.
 *  - maxIter : nombre maximal de vecteurs de Krylov (typiquement ≤ n).
 *  - tol : tolérance pour ||res||/||b||.
 *
 * Retourne le nombre d'itérations effectuées, ou -1 en cas de non-convergence.
 */
int gmres(const double *A, const double *b, double *x, int n,
          int maxIter, double tol)
{
    // Allocation de V : tableau de (maxIter+1) vecteurs de dimension n
    double **V = (double**) malloc((maxIter+1) * sizeof(double*));
    for (int i = 0; i < maxIter+1; i++) {
        V[i] = (double*) calloc(n, sizeof(double));
    }

    // Allocation de H : matrice de taille (maxIter+1) x (maxIter)
    double **H = (double**) malloc((maxIter+1) * sizeof(double*));
    for (int i = 0; i < maxIter+1; i++) {
        H[i] = (double*) calloc(maxIter, sizeof(double));
    }

    // Vecteur g pour le système des moindres carrés
    double *g = (double*) calloc(maxIter+1, sizeof(double));

    // Calcul de r0 = b - A*x
    double *Ax = (double*) malloc(n * sizeof(double));
    matVecProduct(A, x, Ax, n);
    double *r0 = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        r0[i] = b[i] - Ax[i];
    }

    double beta = vecNorm(r0, n);
    double normb = vecNorm(b, n);
    if (normb < 1e-14) {
        normb = 1.0;  // pour éviter la division par zéro si b est nul
    }

    // Si le résidu initial est déjà assez petit, on sort
    if (beta / normb < tol) {
        free(Ax); free(r0);
        for (int i = 0; i < maxIter+1; i++) { free(V[i]); free(H[i]); }
        free(V); free(H); free(g);
        return 0;
    }

    // Initialisation de V[0] = r0 / beta
    for (int i = 0; i < n; i++) {
        V[0][i] = r0[i] / beta;
    }
    g[0] = beta;  // g = [beta, 0, 0, ...]

    int k;
    for (k = 0; k < maxIter; k++) {
        // 1) Calcul de w = A * V[k]
        double *w = (double*) calloc(n, sizeof(double));
        matVecProduct(A, V[k], w, n);

        // 2) Orthogonalisation (Gram-Schmidt modifié)
        for (int j = 0; j <= k; j++) {
            double hij = dotProduct(w, V[j], n);
            H[j][k] = hij;
            for (int i = 0; i < n; i++) {
                w[i] -= hij * V[j][i];
            }
        }
        double normw = vecNorm(w, n);
        H[k+1][k] = normw;

        // 3) Normalisation pour obtenir V[k+1]
        if (normw < 1e-14) {
            // "Lucky breakdown"
            free(w);
            break;
        }
        for (int i = 0; i < n; i++) {
            V[k+1][i] = w[i] / normw;
        }
        free(w);

        // 4) Application d'une rotation de Givens pour trianguler H
        // Ici, nous ne conservons pas les rotations précédentes (simplification)
        // On calcule uniquement la rotation permettant d'annuler H[k+1][k]
        {
            double h1 = H[k][k];
            double h2 = H[k+1][k];
            double denom = sqrt(h1*h1 + h2*h2) + 1e-24;
            double c = h1 / denom;
            double s = -h2 / denom;

            // Appliquer la rotation sur H[k][k] et H[k+1][k]
            H[k][k]   = c * h1 - s * h2;
            H[k+1][k] = 0.0;

            // Mettre à jour g[k] et g[k+1]
            double gk   = g[k];
            double gkp1 = g[k+1];
            g[k]   = c * gk - s * gkp1;
            g[k+1] = s * gk + c * gkp1;
        }

        // Calcul du résidu courant = |g[k+1]|
        double res = fabs(g[k+1]);
        double relRes = res / normb;
        if (relRes < tol) {
            k++;  // on a effectué k+1 itérations
            break;
        }
    }

    // Résolution du système triangulaire H * y = g par substitution arrière
    double *y = (double*) calloc(k+1, sizeof(double));
    for (int i = k-1; i >= 0; i--) {
        double somme = g[i];
        for (int j = i+1; j < k; j++) {
            somme -= H[i][j] * y[j];
        }
        y[i] = somme / H[i][i];
    }

    // Mise à jour de la solution : x = x + V * y
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            x[i] += V[j][i] * y[j];
        }
    }

    // Libération de la mémoire allouée
    free(y);
    free(Ax);
    free(r0);
    for (int i = 0; i < maxIter+1; i++) {
        free(V[i]);
        free(H[i]);
    }
    free(V);
    free(H);
    free(g);

    if (k == maxIter) {
        // Pas de convergence (ou convergence partielle)
        return -1;
    }
    return k; // nombre d'étapes d'Arnoldi utilisées
}

/*-------------------------------------------------------------
 * Fonction main d'exemple
 *-------------------------------------------------------------*/
int main(void) {
    int n = 4;
    double A[16] = {
        1,  3,  2,  0,
        0,  2, -1,  1,
        4,  1,  1,  1,
        0,  1,  2,  3
    };
    double b[4] = { 8, 1, 10, 7 };

    // Initialisation de x à zéro (premier vecteur de départ)
    double x[4] = {0.0, 0.0, 0.0, 0.0};

    double tol = 1e-8;
    int maxIter = 4;  // Pour le GMRES complet, maxIter peut être jusqu'à n

    int iters = gmres(A, b, x, n, maxIter, tol);
    if (iters < 0) {
        printf("GMRES n'a pas convergé en %d itérations.\n", maxIter);
    } else {
        printf("GMRES a convergé en %d itérations.\n", iters);
        printf("Solution x = [");
        for (int i = 0; i < n; i++) {
            printf(" %g", x[i]);
        }
        printf(" ]\n");
    }
    return 0;
}

