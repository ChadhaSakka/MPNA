#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "sparse_matvec.h"

void test_serial(CSRMatrix *csr) {
    double *x = malloc(csr->n * sizeof(double));
    double *y = malloc(csr->n * sizeof(double));
    for (int i = 0; i < csr->n; i++) {
        x[i] = 1.0; // Test vector of all 1s
    }

    double start = MPI_Wtime();
    csr_matvec_serial(csr, x, y);
    double end = MPI_Wtime();

    printf("[Serial] Time: %f seconds\n", end - start);
    free(x);
    free(y);
}

void test_distributed(CSRMatrix *csr, int rank, int size) {
    double *x = malloc(csr->n * sizeof(double));
    double *y = calloc(csr->n, sizeof(double)); // Initialize to 0
    if (rank == 0) {
        for (int i = 0; i < csr->n; i++) {
            x[i] = 1.0; // Test vector of all 1s
        }
    }
    MPI_Bcast(x, csr->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double start = MPI_Wtime();
    csr_matvec_distributed(csr, x, y, rank, size);
    double end = MPI_Wtime();

    if (rank == 0) {
        printf("[Distributed] Time: %f seconds with %d processes\n", end - start, size);
    }
    free(x);
    free(y);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <matrix.mtx>\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    CSRMatrix csr;
    if (rank == 0) {
        read_matrix_market(argv[1], &csr);
        printf("Matrix loaded: n = %d, nnz = %d\n", csr.n, csr.nnz);
    }

    // Broadcast matrix metadata
    MPI_Bcast(&csr.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&csr.nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate and broadcast CSR arrays
    if (rank != 0) {
        csr.values = malloc(csr.nnz * sizeof(double));
        csr.col_ind = malloc(csr.nnz * sizeof(int));
        csr.row_ptr = malloc((csr.n + 1) * sizeof(int));
    }
    MPI_Bcast(csr.row_ptr, csr.n + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(csr.col_ind, csr.nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(csr.values, csr.nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Test serial implementation (only rank 0)
    if (rank == 0 && size == 1) {
        test_serial(&csr);
    }

    // Test distributed implementation
    test_distributed(&csr, rank, size);

    // Test power iteration
    if (rank == 0) {
        printf("Starting power iteration...\n");
    }
    double lambda = power_iteration(&csr, 1000, 1e-6, rank, size);
    if (rank == 0) {
        printf("Dominant eigenvalue: %lf\n", lambda);
    }

    free_csr(&csr);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
