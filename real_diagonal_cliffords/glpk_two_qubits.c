#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <math.h>
#include <time.h>

void Constructor2(float Clifford[32][4],float Final[4][64],int  num_rows, int num_clifford_gates);

int main(void)
{


  glp_prob *lp;

  FILE *fp;

  float (*Clifford)[4]    = malloc(32 * 4 * sizeof *Clifford);
  int i,j;
  float x1,x2,x3,x4;
  int num_qubits = 2;
  double z;


  /*Load Clifford matrices. */
    i=0;j=0;
    fp = fopen("real_diagonal_cliffords/two_qubit_real_diagonal_clifford_groups.txt","r");
    while(fscanf(fp ,"%f\t%f\t%f\t%f", &x1, &x2, &x3, &x4)!=EOF)
      {
        Clifford[i][0]=x1;
        Clifford[i][1]=x2;
        Clifford[i][2]=x3;
        Clifford[i][3]=x4;
        i=i+1;
      }

    fclose(fp);
    i=0;j=0;

  /*Define Variables to solve the problem*/
  int num_clifford_gates= 32;

  


  int num_rows =  pow(2, num_qubits);
  int num_cols = 2*num_clifford_gates;



  //float (*Final)[num_cols]   = malloc(num_rows * num_cols * sizeof *Final);
  float (*Final)[num_cols] = malloc(num_rows * sizeof(*Final));

  int    *row    = malloc((num_rows*num_cols+1) * sizeof *row);
  int    *column = malloc((num_rows*num_cols+1) * sizeof *column);
  double *values = malloc((num_rows*num_cols+1) * sizeof *values);


  /* Obtain the coefficient matrix*/
    Constructor2(Clifford ,Final, num_rows, num_clifford_gates);




  lp = glp_create_prob();
  glp_set_prob_name(lp, "CZ");
  glp_set_obj_dir(lp, GLP_MIN);


  glp_add_rows(lp, num_rows); /*rows, elements of the matrix to be decomposed (flattened)*/

   /*Real --------------------------------*/
   for (i = 1; i < num_rows + 1; i++){
    if (i == num_rows){
      glp_set_row_bnds(lp, i, GLP_FX, -1.0, -1.0);
    } else{
      glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0);
    }
   }



  glp_add_cols(lp, (2*num_clifford_gates));/*Número de columnas que se añaden. Vamos, variables libres que tiene el problema*/
  /*Dimension=#elementosgrupo*2(real e imaginaria)*/
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


void Constructor2(float Clifford[32][4],float Final[4][64],int  num_rows, int num_clifford_gates)
{
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

