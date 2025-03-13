#ifndef SOLVEUR_H
#define SOLVEUR_H

/**
 * @brief Résout un système linéaire Ax = b en utilisant l'interface IJ de Hypre.
 *
 * Les matrices sont supposées stockées au format CSR.
 *
 * @param N Taille du système (nombre d'inconnues).
 * @param A_values Tableau contenant les valeurs non nulles de la matrice.
 * @param A_row_indices Tableau des indices de début de chaque ligne.
 * @param A_col_indices Tableau des indices de colonnes associés aux valeurs.
 * @param nnz Nombre de valeurs non nulles.
 * @param b Vecteur second membre.
 * @param x Vecteur solution (sortie).
 * @return int 0 en cas de succès, code d'erreur sinon.
 */
int solve_linear_system(int N, double *A_values, int *A_row_indices, int *A_col_indices, int nnz,
                          double *b, double *x);

#endif // SOLVEUR_H

