#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <math.h>
#include <time.h>

/*
 2 qubits: num_clifford_gates = 32
 3 qubits: num_clifford_gates = 64
 4 qubits: num_clifford_gates = 1024
 5 qubits: num_clifford_gates = 32768
 6 qubits: num-clifford_gates = 2097152

*/
int num_clifford_gates= 2097152;
int num_qubits = 6;
const char *file_path = "real_diagonal_cliffords/six_qubit_real_diagonal_clifford_groups.txt";



void CoefficientMatrix(float Clifford[num_clifford_gates][1 << num_qubits], float Final[1 << num_qubits][2 * num_clifford_gates]);

int main(void)
{

  int num_cols = 2 * num_clifford_gates;
  int num_rows = 1 << num_qubits;
  double z;


  glp_prob *lp;

  FILE *fp;

  float (*Clifford)[1 << num_qubits] = malloc(num_clifford_gates * (1 << num_qubits) * sizeof *Clifford);
  int i,j;
  
  /*Read the Clifford matrix from file*/

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

  /*Define Variables to solve the problem*/
  //float (*Final)[num_cols]   = malloc(num_rows * num_cols * sizeof *Final);
  float (*Final)[num_cols] = malloc(num_rows * sizeof(*Final));

  int    *row    = malloc((num_rows*num_cols+1) * sizeof *row);
  int    *column = malloc((num_rows*num_cols+1) * sizeof *column);
  double *values = malloc((num_rows*num_cols+1) * sizeof *values);


  /* Obtain the coefficient matrix*/
  CoefficientMatrix(Clifford ,Final);




  lp = glp_create_prob();
  glp_set_prob_name(lp, "Multi CZ gate");
  glp_set_obj_dir(lp, GLP_MIN);


  glp_add_rows(lp, num_rows); /*rows, elements of the matrix to be decomposed (flattened)*/

   /*rows are the constraints. Here, it is the multi CZ gate to be decomposed. Ajust this for the real diagonal gate you want to decompose*/
   for (i = 1; i < num_rows + 1; i++){
    if (i == num_rows){
      glp_set_row_bnds(lp, i, GLP_FX, -1.0, -1.0);
    } else{
      glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0);
    }
   }



  glp_add_cols(lp, (2*num_clifford_gates));/*Columns are the variables of the problem. */

  for ( i = 1; i < (2*num_clifford_gates + 1); i++)
  {
    glp_set_col_bnds(lp,(i),GLP_LO,0.0,0.0);
    glp_set_obj_coef(lp, (i), 1.0);
  }
  i=0;

  int k,l;

  k=1;
  l=1;
  for (j = 1; j < num_rows+1; j++)
  {
    for ( i = 1; i < num_cols+1; i++)
    {


      row[l]=j;
      column[l]=i;
      values[l]=Final[j-1][i-1];

      l=l+1;



    }

  }

  clock_t start, end;
  double cpu_time_used;

  glp_load_matrix(lp, num_cols*num_rows, row, column, values);
  start = clock();
  glp_simplex(lp, NULL);
  end = clock();

  z = glp_get_obj_val(lp);
  double z_squared = z*z;
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("The optimal z value is: %.10f\n", z);
  printf("The optimal z value squared is: %.10f\n", z_squared);
  printf("Time taken: %f seconds\n", cpu_time_used);


  glp_delete_prob(lp);
  glp_free_env();
  return 0;
  }


void CoefficientMatrix(float Clifford[num_clifford_gates][1 << num_qubits], float Final[1 << num_qubits][2 * num_clifford_gates])
{
  int num_rows = 1 << num_qubits;
  int j, i;

  for (j = 0; j < num_rows; j++)
  {
    for (i = 0; i < num_clifford_gates; i++)
    {
      Final[j][i] = Clifford[i][j];
      Final[j][i + num_clifford_gates] = -Clifford[i][j];
    }
  }
}
