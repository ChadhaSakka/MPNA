#ifndef SPARSE_MATVEC_H
#define SPARSE_MATVEC_H

#include <mpi.h>

// CSR Matrix structure
typedef struct {
    int n;        // Matrix dimension (n x n)
    int nnz;      // Number of non-zero elements
    int *row_ptr; // Row pointers
    int *col_ind; // Column indices
    double *values; // Non-zero values
} CSRMatrix;

// Function declarations
void read_matrix_market(const char *filename, CSRMatrix *csr);
void csr_matvec_serial(CSRMatrix *csr, double *x, double *y);
void csr_matvec_distributed(CSRMatrix *csr, double *x, double *y, int rank, int size);
void normalize(double *x, int n);
double power_iteration(CSRMatrix *csr, int max_iter, double tol, int rank, int size);
void free_csr(CSRMatrix *csr);

#endif // SPARSE_MATVEC_H
