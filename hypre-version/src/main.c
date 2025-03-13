#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"
#include "diffusion.h"
#include "utils.h"

// Méthode explicite
void solve_explicit(double *u, int N, const parameters_t *params, double dx, double gamma, int *iter_out, double *final_diff) {
    double *u_new = alloc_vector(N + 1);
    double tol = 1e-6;
    int max_iter = 10000;
    int iter = 0;
    double diff = tol + 1.0;

    while (diff > tol && iter < max_iter) {
        double dt = compute_dt(u, N, params, dx, gamma);
        for (int i = 0; i <= N; i++) {
            double kappa_im = (i > 0) ? 0.5 * (compute_kappa(u[i-1], params) + compute_kappa(u[i], params)) : 0.0;
            double kappa_ip = (i < N) ? 0.5 * (compute_kappa(u[i], params) + compute_kappa(u[i+1], params)) : 0.0;
            double diffusion = (i > 0 && i < N) ? (kappa_ip * (u[i+1] - u[i]) - kappa_im * (u[i] - u[i-1])) / (dx * dx) : 0.0;
            if (i == 0) diffusion = (kappa_ip * (u[1] - u[0])) / (dx * dx); // Neumann
            if (i == N) diffusion = 0.0; // Dirichlet imposé après
            double x = i * dx;
            u_new[i] = u[i] + dt * (diffusion - params->sigma * (pow(u[i], 4) - 1.0) + source_Q(x, params));
        }
        u_new[N] = 1.0; // Dirichlet
        diff = 0.0;
        for (int i = 0; i <= N; i++) {
            diff += fabs(u_new[i] - u[i]);
            u[i] = u_new[i];
        }
        diff /= (N + 1);
        iter++;
    }
    *iter_out = iter;
    *final_diff = diff;
    free_vector(u_new);
}

// Méthode implicite avec HYPRE
void solve_implicit(double *u, int N, const parameters_t *params, double dx, double gamma, int *iter_out, double *final_diff) {
    double tol = 1e-6;
    int max_iter = 1000;
    int iter = 0;
    double diff = tol + 1.0;
    double *u_new = alloc_vector(N + 1);
    double *Q = alloc_vector(N + 1);
    for (int i = 0; i <= N; i++) {
        Q[i] = source_Q(i * dx, params);
    }

    HYPRE_IJMatrix A;
    HYPRE_IJVector b, x;
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b, par_x;
    HYPRE_Solver solver, precond;

    while (diff > tol && iter < max_iter) {
        double dt = compute_dt(u, N, params, dx, gamma);

        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N, 0, N, &A);
        HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N, &b);
        HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N, &x);
        HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(x);

        for (int i = 0; i <= N; i++) {
            double kappa_im = (i > 0) ? 0.5 * (compute_kappa(u[i-1], params) + compute_kappa(u[i], params)) : 0.0;
            double kappa_ip = (i < N) ? 0.5 * (compute_kappa(u[i], params) + compute_kappa(u[i+1], params)) : 0.0;
            if (i == 0) { // Neumann
                double a = 1.0 + dt * (kappa_ip / (dx * dx));
                double b_val = -dt * (kappa_ip / (dx * dx));
                int cols[2] = {i, i+1};
                double values[2] = {a, b_val};
                HYPRE_IJMatrixSetValues(A, 1, &(int){2}, &i, cols, values);
                double rhs = u[i] + dt * (Q[i] + params->sigma);
                HYPRE_IJVectorSetValues(b, 1, &i, &rhs);
            } else if (i == N) { // Dirichlet
                int cols[1] = {i};
                double values[1] = {1.0};
                HYPRE_IJMatrixSetValues(A, 1, &(int){1}, &i, cols, values);
                double rhs = 1.0;
                HYPRE_IJVectorSetValues(b, 1, &i, &rhs);
            } else { // Points internes
                double a = 1.0 + dt * ((kappa_im + kappa_ip) / (dx * dx) + params->sigma * pow(u[i], 3));
                double b_val = -dt * (kappa_ip / (dx * dx));
                double c_val = -dt * (kappa_im / (dx * dx));
                int cols[3] = {i-1, i, i+1};
                double values[3] = {c_val, a, b_val};
                HYPRE_IJMatrixSetValues(A, 1, &(int){3}, &i, cols, values);
                double rhs = u[i] + dt * (Q[i] + params->sigma);
                HYPRE_IJVectorSetValues(b, 1, &i, &rhs);
            }
        }

        HYPRE_IJMatrixAssemble(A);
        HYPRE_IJMatrixGetObject(A, (void**)&par_A);
        HYPRE_IJVectorAssemble(b);
        HYPRE_IJVectorGetObject(b, (void**)&par_b);
        HYPRE_IJVectorAssemble(x);
        HYPRE_IJVectorGetObject(x, (void**)&par_x);

        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetTol(precond, 1e-8);
        HYPRE_BoomerAMGSetMaxIter(precond, 100);

        HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_GMRESSetMaxIter(solver, 100);
        HYPRE_GMRESSetTol(solver, 1e-8);
        HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);
        HYPRE_GMRESSetup(solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);
        HYPRE_GMRESSolve(solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);

        for (int i = 0; i <= N; i++) {
            HYPRE_IJVectorGetValues(x, 1, &i, &u_new[i]);
        }

        diff = 0.0;
        for (int i = 0; i <= N; i++) {
            diff += fabs(u_new[i] - u[i]);
            u[i] = u_new[i];
        }
        diff /= (N + 1);
        iter++;

        HYPRE_IJMatrixDestroy(A);
        HYPRE_IJVectorDestroy(b);
        HYPRE_IJVectorDestroy(x);
        HYPRE_ParCSRGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    *iter_out = iter;
    *final_diff = diff;
    free_vector(u_new);
    free_vector(Q);
}

// Méthode de Newton-Raphson avec HYPRE
void solve_newton(double *u, int N, const parameters_t *params, double dx, int *iter_out, double *final_diff) {
    int max_iter = 300;
    double tol = 1e-6;
    int iter = 0;
    double *F = alloc_vector(N + 1);
    double *u_new = alloc_vector(N + 1);

    HYPRE_IJMatrix J;
    HYPRE_IJVector b, delta;
    HYPRE_ParCSRMatrix par_J;
    HYPRE_ParVector par_b, par_delta;
    HYPRE_Solver solver, precond;

    while (iter < max_iter) {
        for (int i = 0; i <= N; i++) {
            double kappa_im = (i > 0) ? 0.5 * (compute_kappa(u[i-1], params) + compute_kappa(u[i], params)) : 0.0;
            double kappa_ip = (i < N) ? 0.5 * (compute_kappa(u[i], params) + compute_kappa(u[i+1], params)) : 0.0;
            double x = i * dx;
            if (i == 0) {
                F[i] = u[0] - u[1]; // Neumann: du/dx = 0
            } else if (i == N) {
                F[i] = u[N] - 1.0; // Dirichlet: u = 1
            } else {
                F[i] = -(kappa_ip * (u[i+1] - u[i]) - kappa_im * (u[i] - u[i-1])) / (dx * dx)
                       + params->sigma * (pow(u[i], 4) - 1.0) - source_Q(x, params);
            }
        }

        double res_norm = norm_vector(F, N + 1);
        if (res_norm < tol) break;

        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N, 0, N, &J);
        HYPRE_IJMatrixSetObjectType(J, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(J);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N, &b);
        HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b);
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N, &delta);
        HYPRE_IJVectorSetObjectType(delta, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(delta);

        for (int i = 0; i <= N; i++) {
            double kappa_im = (i > 0) ? 0.5 * (compute_kappa(u[i-1], params) + compute_kappa(u[i], params)) : 0.0;
            double kappa_ip = (i < N) ? 0.5 * (compute_kappa(u[i], params) + compute_kappa(u[i+1], params)) : 0.0;
            if (i == 0) {
                int cols[2] = {i, i+1};
                double values[2] = {1.0, -1.0};
                HYPRE_IJMatrixSetValues(J, 1, &(int){2}, &i, cols, values);
                double temp = -F[i];
                HYPRE_IJVectorSetValues(b, 1, &i, &temp);
            } else if (i == N) {
                int cols[1] = {i};
                double values[1] = {1.0};
                HYPRE_IJMatrixSetValues(J, 1, &(int){1}, &i, cols, values);
                double temp = -F[i];
                HYPRE_IJVectorSetValues(b, 1, &i, &temp);
            } else {
                double a = (kappa_im + kappa_ip) / (dx * dx) + 4.0 * params->sigma * pow(u[i], 3);
                double b_val = -kappa_ip / (dx * dx);
                double c_val = -kappa_im / (dx * dx);
                int cols[3] = {i-1, i, i+1};
                double values[3] = {c_val, a, b_val};
                HYPRE_IJMatrixSetValues(J, 1, &(int){3}, &i, cols, values);
                double temp = -F[i];
                HYPRE_IJVectorSetValues(b, 1, &i, &temp);
            }
        }

        HYPRE_IJMatrixAssemble(J);
        HYPRE_IJMatrixGetObject(J, (void**)&par_J);
        HYPRE_IJVectorAssemble(b);
        HYPRE_IJVectorGetObject(b, (void**)&par_b);
        HYPRE_IJVectorAssemble(delta);
        HYPRE_IJVectorGetObject(delta, (void**)&par_delta);

        HYPRE_BoomerAMGCreate(&precond);
        HYPRE_BoomerAMGSetTol(precond, 1e-8);
        HYPRE_BoomerAMGSetMaxIter(precond, 100);

        HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_GMRESSetMaxIter(solver, 100);
        HYPRE_GMRESSetTol(solver, 1e-8);
        HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                              (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);
        HYPRE_GMRESSetup(solver, (HYPRE_Matrix)par_J, (HYPRE_Vector)par_b, (HYPRE_Vector)par_delta);
        HYPRE_GMRESSolve(solver, (HYPRE_Matrix)par_J, (HYPRE_Vector)par_b, (HYPRE_Vector)par_delta);

        for (int i = 0; i <= N; i++) {
            double delta_val;
            HYPRE_IJVectorGetValues(delta, 1, &i, &delta_val);
            u_new[i] = u[i] + delta_val;
            if (u_new[i] < 1e-10) u_new[i] = 1e-10; // Éviter les valeurs négatives
            u[i] = u_new[i];
        }

        iter++;
        HYPRE_IJMatrixDestroy(J);
        HYPRE_IJVectorDestroy(b);
        HYPRE_IJVectorDestroy(delta);
        HYPRE_ParCSRGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    *iter_out = iter;
    *final_diff = norm_vector(F, N + 1);
    free_vector(F);
    free_vector(u_new);
}
 // Ajoutez cette fonction pour sauvegarder les solutions dans un fichier CSV
void write_solutions_to_file(const char *filename, double *u_explicit, double *u_implicit, double *u_newton, int N, double dx) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Erreur d'ouverture du fichier %s\n", filename);
        return;
    }
    // Entête du fichier CSV
    fprintf(fp, "x,u_explicit,u_implicit,u_newton\n");
    for (int i = 0; i <= N; i++) {
        double x = i * dx;
        fprintf(fp, "%f,%f,%f,%f\n", x, u_explicit[i], u_implicit[i], u_newton[i]);
    }
    fclose(fp);
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    HYPRE_Init();

    parameters_t params = {0.01, 0.1, 1.0, 0.5}; // kappa0, sigma, beta, q
    int N = 100;
    double dx = 1.0 / N;
    double gamma = 0.1;
    int iter_explicit, iter_implicit, iter_newton;
    double diff_explicit, diff_implicit, diff_newton;
    double t_explicit, t_implicit, t_newton;

    double *u_explicit = alloc_vector(N + 1);
    double *u_implicit = alloc_vector(N + 1);
    double *u_newton = alloc_vector(N + 1);
    discretize_diffusion(u_explicit, N, &params, dx);
    discretize_diffusion(u_implicit, N, &params, dx);
    discretize_diffusion(u_newton, N, &params, dx);

    printf("=== Résolution de l'équation de diffusion non linéaire ===\n");
    printf("Paramètres: κ₀ = %g, σ = %g, β = %g, q = %g\n", params.kappa0, params.sigma, params.beta, params.q);
    printf("Grille: N = %d, dx = %g, γ = %g\n\n", N, dx, gamma);

    double start = MPI_Wtime();
    solve_explicit(u_explicit, N, &params, dx, gamma, &iter_explicit, &diff_explicit);
    t_explicit = MPI_Wtime() - start;
    printf("Méthode explicite: %d itérations, erreur = %e, temps = %g s\n", iter_explicit, diff_explicit, t_explicit);

    start = MPI_Wtime();
    solve_implicit(u_implicit, N, &params, dx, gamma, &iter_implicit, &diff_implicit);
    t_implicit = MPI_Wtime() - start;
    printf("Méthode implicite: %d itérations, erreur = %e, temps = %g s\n", iter_implicit, diff_implicit, t_implicit);

    start = MPI_Wtime();
    solve_newton(u_newton, N, &params, dx, &iter_newton, &diff_newton);
    t_newton = MPI_Wtime() - start;
    printf("Méthode de Newton: %d itérations, erreur = %e, temps = %g s\n", iter_newton, diff_newton, t_newton);

    printf("\nSolution finale (quelques points) - Méthode explicite:\n");
    for (int i = 0; i <= N; i += N/5) printf("x = %.4f, u = %.6f\n", i * dx, u_explicit[i]);
    printf("\nSolution finale - Méthode implicite:\n");
    for (int i = 0; i <= N; i += N/5) printf("x = %.4f, u = %.6f\n", i * dx, u_implicit[i]);
    printf("\nSolution finale - Méthode de Newton:\n");
    for (int i = 0; i <= N; i += N/5) printf("x = %.4f, u = %.6f\n", i * dx, u_newton[i]);

    // Sauvegarde des solutions dans un fichier pour post-traitement
    write_solutions_to_file("solutions.csv", u_explicit, u_implicit, u_newton, N, dx);
  
    free_vector(u_explicit);
    free_vector(u_implicit);
    free_vector(u_newton);

    HYPRE_Finalize();
    MPI_Finalize();
    return 0;
}
