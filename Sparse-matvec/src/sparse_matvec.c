#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sparse_matvec.h"

#define MAX_LINE 1024

void read_matrix_market(const char *filename, CSRMatrix *csr) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char line[MAX_LINE];
    while (fgets(line, MAX_LINE, file) && line[0] == '%'); // Skip comments
    int rows, cols, nnz;
    sscanf(line, "%d %d %d", &rows, &cols, &nnz);

    csr->n = rows;
    csr->nnz = nnz;
    csr->values = malloc(nnz * sizeof(double));
    csr->col_ind = malloc(nnz * sizeof(int));
    csr->row_ptr = calloc(rows + 1, sizeof(int));

    int *row_counts = calloc(rows, sizeof(int));
    int row, col;
    double value;
    for (int i = 0; i < nnz; i++) {
        fscanf(file, "%d %d %lf", &row, &col, &value);
        row--; col--; // 0-based indexing
        row_counts[row]++;
    }

    for (int i = 1; i <= rows; i++) {
        csr->row_ptr[i] = csr->row_ptr[i - 1] + row_counts[i - 1];
    }

    rewind(file);
    while (fgets(line, MAX_LINE, file) && line[0] == '%'); // Skip comments again
    int *offsets = calloc(rows, sizeof(int));
    for (int i = 0; i < nnz; i++) {
        fscanf(file, "%d %d %lf", &row, &col, &value);
        row--; col--;
        int idx = csr->row_ptr[row] + offsets[row];
        csr->col_ind[idx] = col;
        csr->values[idx] = value;
        offsets[row]++;
    }

    fclose(file);
    free(row_counts);
    free(offsets);
}

void csr_matvec_serial(CSRMatrix *csr, double *x, double *y) {
    for (int i = 0; i < csr->n; i++) {
        y[i] = 0.0;
        for (int j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++) {
            y[i] += csr->values[j] * x[csr->col_ind[j]];
        }
    }
}

void csr_matvec_distributed(CSRMatrix *csr, double *x, double *y, int rank, int size) {
    int block_size = csr->n / size;
    int start_row = rank * block_size;
    int end_row = (rank == size - 1) ? csr->n : start_row + block_size;

    // Local computation
    for (int i = start_row; i < end_row; i++) {
        y[i] = 0.0;
        for (int j = csr->row_ptr[i]; j < csr->row_ptr[i + 1]; j++) {
            y[i] += csr->values[j] * x[csr->col_ind[j]];
        }
    }

    // Gather results to all processes
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DOUBLE, y, block_size, MPI_DOUBLE, MPI_COMM_WORLD);
}

void normalize(double *x, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += x[i] * x[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) {
        x[i] /= norm;
    }
}

double power_iteration(CSRMatrix *csr, int max_iter, double tol, int rank, int size) {
    double *x = malloc(csr->n * sizeof(double));
    double *y = malloc(csr->n * sizeof(double));

    // Initialize x with uniform values (only rank 0 initializes, then broadcasts)
    if (rank == 0) {
        for (int i = 0; i < csr->n; i++) {
            x[i] = 1.0 / csr->n;
        }
        normalize(x, csr->n);
    }
    MPI_Bcast(x, csr->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double lambda = 0.0, lambda_old = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        csr_matvec_distributed(csr, x, y, rank, size);
        if (rank == 0) {
            normalize(y, csr->n);
            lambda_old = lambda;
            lambda = 0.0;
            for (int i = 0; i < csr->n; i++) {
                lambda += x[i] * y[i];
            }
            if (iter % 100 == 0) {
                printf("Iteration %d: lambda = %lf\n", iter, lambda);
            }
            if (fabs(lambda - lambda_old) < tol) {
                printf("Converged at iteration %d with lambda = %lf\n", iter, lambda);
                break;
            }
            memcpy(x, y, csr->n * sizeof(double));
        }
        MPI_Bcast(x, csr->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    free(x);
    free(y);
    return lambda;
}

void free_csr(CSRMatrix *csr) {
    free(csr->values);
    free(csr->col_ind);
    free(csr->row_ptr);
}
