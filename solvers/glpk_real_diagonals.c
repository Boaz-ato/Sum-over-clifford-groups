/*
 * glpk_real_diagonals.c
 *
 * Solves the Clifford decomposition linear programme using the GLPK simplex solver.
 *
 * Problem:
 *   Given a target diagonal gate b ∈ {±1}^(2^n), find non-negative coefficients x
 *   that minimise sum(x) subject to A*x = b, where the columns of A are the diagonals
 *   of all real diagonal Clifford group elements (and their negatives).
 *
 * To configure the qubit count and target gate, edit the global variables at the top
 * of this file before compiling:
 *
 *   num_clifford_gates:
 *     2 qubits:  32
 *     3 qubits:  64
 *     4 qubits:  1024
 *     5 qubits:  32768
 *     6 qubits:  2097152
 *
 *   num_qubits: set to match the chosen num_clifford_gates above.
 *   file_path:  path to the corresponding dataset file.
 *
 * To select a different target gate, edit the row bounds loop in main() to match
 * the desired diagonal vector.
 *
 * Compile (from the project root):
 *   gcc -O2 -o glpk_solver solvers/glpk_real_diagonals.c -lglpk -lm
 *
 * Run:
 *   ./glpk_solver
 */

#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <math.h>
#include <time.h>

int num_clifford_gates = 2097152;
int num_qubits = 6;
const char *file_path = "datasets/real_diagonals/six_qubit_real_diagonal_clifford_groups.txt";

/*
 * CoefficientMatrix - Build the LP coefficient matrix from the Clifford group.
 *
 * Each Clifford group element contributes two columns: the original diagonal
 * and its negation, so that non-negative LP variables can represent both
 * positive and negative coefficients.
 *
 * Args:
 *   Clifford: Input array of shape [num_clifford_gates][2^num_qubits] containing
 *             the diagonal entries of each group element.
 *   Final:    Output array of shape [2^num_qubits][2*num_clifford_gates] to be
 *             filled with the coefficient matrix.
 */
void CoefficientMatrix(float Clifford[num_clifford_gates][1 << num_qubits],
                       float Final[1 << num_qubits][2 * num_clifford_gates]);

int main(void)
{
    int num_cols = 2 * num_clifford_gates;
    int num_rows = 1 << num_qubits;
    double z;

    glp_prob *lp;
    FILE *fp;

    float (*Clifford)[1 << num_qubits] = malloc(num_clifford_gates * (1 << num_qubits) * sizeof *Clifford);
    int i, j;

    /* Read the Clifford group diagonal entries from file */
    float *x_values = malloc((1 << num_qubits) * sizeof(float));
    if (x_values == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }

    fp = fopen(file_path, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open file %s\n", file_path);
        free(x_values);
        return 1;
    }

    i = 0;
    while (1) {
        for (j = 0; j < (1 << num_qubits); j++) {
            if (j < (1 << num_qubits) - 1) {
                if (fscanf(fp, "%f\t", &x_values[j]) != 1) {
                    goto end_read;
                }
            } else {
                if (fscanf(fp, "%f", &x_values[j]) != 1) {
                    goto end_read;
                }
            }
        }

        for (j = 0; j < (1 << num_qubits); j++) {
            Clifford[i][j] = x_values[j];
        }
        i++;
    }

end_read:
    fclose(fp);
    free(x_values);

    /* Allocate LP coefficient matrix and index arrays */
    float (*Final)[num_cols] = malloc(num_rows * sizeof(*Final));
    int    *row    = malloc((num_rows * num_cols + 1) * sizeof *row);
    int    *column = malloc((num_rows * num_cols + 1) * sizeof *column);
    double *values = malloc((num_rows * num_cols + 1) * sizeof *values);

    CoefficientMatrix(Clifford, Final);

    /* Set up the GLPK problem */
    lp = glp_create_prob();
    glp_set_prob_name(lp, "Multi CZ gate");
    glp_set_obj_dir(lp, GLP_MIN);

    /* Add rows (constraints). Each row corresponds to one diagonal entry of the
     * target gate. Adjust bounds here to match a different target gate:
     *   - All entries +1 except the last (-1) gives a C^(n-1)Z gate.
     *   - Other patterns correspond to tensor products or circuits of gates. */
    glp_add_rows(lp, num_rows);
    for (i = 1; i < num_rows + 1; i++) {
        if (i == num_rows) {
            glp_set_row_bnds(lp, i, GLP_FX, -1.0, -1.0);
        } else {
            glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0);
        }
    }

    /* Add columns (decision variables), one per Clifford element (positive and
     * negative), with non-negativity constraints and unit cost coefficient. */
    glp_add_cols(lp, 2 * num_clifford_gates);
    for (i = 1; i < (2 * num_clifford_gates + 1); i++) {
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, i, 1.0);
    }

    /* Fill the sparse matrix representation */
    int k = 1;
    for (j = 1; j < num_rows + 1; j++) {
        for (i = 1; i < num_cols + 1; i++) {
            row[k]    = j;
            column[k] = i;
            values[k] = Final[j - 1][i - 1];
            k++;
        }
    }

    /* Load the matrix and solve */
    clock_t start, end;
    double cpu_time_used;

    glp_load_matrix(lp, num_cols * num_rows, row, column, values);
    start = clock();
    glp_simplex(lp, NULL);
    end = clock();

    z = glp_get_obj_val(lp);
    double z_squared = z * z;
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("The optimal z value is: %.10f\n", z);
    printf("The optimal z value squared is: %.10f\n", z_squared);
    printf("Time taken: %f seconds\n", cpu_time_used);

    /* Report non-zero Clifford group elements in the solution */
    int nonzero_count = 0;
    printf("\nNonzero Clifford group elements:\n");
    for (i = 1; i <= num_cols; i++) {
        double val = glp_get_col_prim(lp, i);
        if (fabs(val) > 1e-10) {
            nonzero_count++;
            printf("x%d = %.10f\n", i, val);
        }
    }
    printf("\nTotal nonzero Clifford group elements: %d\n", nonzero_count);

    glp_delete_prob(lp);
    glp_free_env();
    return 0;
}


void CoefficientMatrix(float Clifford[num_clifford_gates][1 << num_qubits],
                       float Final[1 << num_qubits][2 * num_clifford_gates])
{
    int num_rows = 1 << num_qubits;
    int j, i;

    for (j = 0; j < num_rows; j++) {
        for (i = 0; i < num_clifford_gates; i++) {
            Final[j][i]                     =  Clifford[i][j];
            Final[j][i + num_clifford_gates] = -Clifford[i][j];
        }
    }
}
