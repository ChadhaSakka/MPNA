#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>           // Ajoutez cette inclusion pour MPI
#include "solveur.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    // Initialiser MPI
    MPI_Init(&argc, &argv);

    // Test simple : résolution d'un système linéaire identité de taille N
    int N = 10;
    double *b = alloc_vector(N);
    double *x = alloc_vector(N);
    double *A_values = alloc_vector(N); // Une valeur par ligne pour la diagonale
    int *A_row_indices = (int *) malloc(N * sizeof(int));
    int *A_col_indices = (int *) malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) {
        b[i] = 1.0;       // Second membre constant
        A_values[i] = 1.0;  // Diagonale = 1
        A_row_indices[i] = i;
        A_col_indices[i] = i;
    }
    int nnz = N; // Nombre de coefficients non nuls

    int ierr = solve_linear_system(N, A_values, A_row_indices, A_col_indices, nnz, b, x);
    if (ierr == 0) {
        printf("Solution du système linéaire (matrice identité):\n");
        print_vector(x, N);
    } else {
        printf("Erreur dans le solveur linéaire, code: %d\n", ierr);
    }

    free_vector(b);
    free_vector(x);
    free(A_values);
    free(A_row_indices);
    free(A_col_indices);

    // Finaliser MPI
    MPI_Finalize();

    return 0;
}

