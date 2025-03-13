#include <stdio.h>
#include <stdlib.h>
#include "solveur.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"

/**
 * @brief Résout un système linéaire en utilisant l'interface IJ de Hypre.
 *
 * Ici, la matrice est construite comme une matrice identité pour illustrer l'utilisation.
 * Dans une application réelle, il faudra construire la matrice à partir des coefficients
 * de la discrétisation.
 *
 * @param N Nombre d'inconnues.
 * @param A_values Valeurs non utilisées dans cet exemple.
 * @param A_row_indices Indices de lignes non utilisés.
 * @param A_col_indices Indices de colonnes non utilisés.
 * @param nnz Nombre d'entrées non nulles.
 * @param b Vecteur du second membre.
 * @param x Vecteur solution (résultat).
 * @return int Code d'erreur (0 si succès).
 */
int solve_linear_system(int N, double *A_values, int *A_row_indices, int *A_col_indices, int nnz,
                          double *b, double *x) {
    int ierr;
    int ilower = 0, iupper = N - 1;
    HYPRE_IJMatrix A;
    HYPRE_IJVector b_vec, x_vec;
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b, par_x;
    
    /* Création et initialisation de la matrice IJ */
    if ((ierr = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A))) {
        fprintf(stderr, "Erreur lors de la création de la matrice IJ\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR))) {
        fprintf(stderr, "Erreur lors de la définition du type de matrice\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJMatrixInitialize(A))) {
        fprintf(stderr, "Erreur lors de l'initialisation de la matrice IJ\n");
        return ierr;
    }

    /* Exemple de remplissage : matrice identité */
    for (int i = 0; i < N; i++) {
        int ncols = 1;
        int col = i;
        double value = 1.0;
        if ((ierr = HYPRE_IJMatrixSetValues(A, 1, &ncols, &i, &col, &value))) {
            fprintf(stderr, "Erreur lors du remplissage de la matrice IJ\n");
            return ierr;
        }
    }
    if ((ierr = HYPRE_IJMatrixAssemble(A))) {
        fprintf(stderr, "Erreur lors de l'assemblage de la matrice IJ\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJMatrixGetObject(A, (void **) &par_A))) {
        fprintf(stderr, "Erreur lors de la récupération de l'objet matrice\n");
        return ierr;
    }

    /* Création et initialisation des vecteurs b et x */
    if ((ierr = HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b_vec))) {
        fprintf(stderr, "Erreur lors de la création du vecteur b\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorSetObjectType(b_vec, HYPRE_PARCSR))) {
        fprintf(stderr, "Erreur lors du paramétrage du vecteur b\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorInitialize(b_vec))) {
        fprintf(stderr, "Erreur lors de l'initialisation du vecteur b\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x_vec))) {
        fprintf(stderr, "Erreur lors de la création du vecteur x\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorSetObjectType(x_vec, HYPRE_PARCSR))) {
        fprintf(stderr, "Erreur lors du paramétrage du vecteur x\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorInitialize(x_vec))) {
        fprintf(stderr, "Erreur lors de l'initialisation du vecteur x\n");
        return ierr;
    }

    /* Remplissage du vecteur b */
    for (int i = 0; i < N; i++) {
        if ((ierr = HYPRE_IJVectorSetValues(b_vec, 1, &i, &b[i]))) {
            fprintf(stderr, "Erreur lors du remplissage du vecteur b\n");
            return ierr;
        }
    }
    if ((ierr = HYPRE_IJVectorAssemble(b_vec))) {
        fprintf(stderr, "Erreur lors de l'assemblage du vecteur b\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorAssemble(x_vec))) {
        fprintf(stderr, "Erreur lors de l'assemblage du vecteur x\n");
        return ierr;
    }
    
    if ((ierr = HYPRE_IJVectorGetObject(b_vec, (void **) &par_b))) {
        fprintf(stderr, "Erreur lors de la récupération du vecteur par b\n");
        return ierr;
    }
    if ((ierr = HYPRE_IJVectorGetObject(x_vec, (void **) &par_x))) {
        fprintf(stderr, "Erreur lors de la récupération du vecteur par x\n");
        return ierr;
    }

    /* Configuration et appel du solveur BoomerAMG */
    HYPRE_Solver solver;
    if ((ierr = HYPRE_BoomerAMGCreate(&solver))) {
        fprintf(stderr, "Erreur lors de la création du solveur BoomerAMG\n");
        return ierr;
    }
    HYPRE_BoomerAMGSetPrintLevel(solver, 1);
    HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.25);
    if ((ierr = HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x))) {
        fprintf(stderr, "Erreur lors du setup du solveur BoomerAMG\n");
        return ierr;
    }
    if ((ierr = HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x))) {
        fprintf(stderr, "Erreur lors de la résolution avec BoomerAMG\n");
        return ierr;
    }

    /* Récupération de la solution */
    for (int i = 0; i < N; i++) {
        if ((ierr = HYPRE_IJVectorGetValues(x_vec, 1, &i, &x[i]))) {
            fprintf(stderr, "Erreur lors de la récupération des valeurs de x\n");
            return ierr;
        }
    }

    /* Libération des objets Hypre */
    HYPRE_BoomerAMGDestroy(solver);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b_vec);
    HYPRE_IJVectorDestroy(x_vec);

    return ierr;
}

