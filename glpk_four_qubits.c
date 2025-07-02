#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <math.h>
#include <time.h>

void Constructor1(float Clifford[524288][16], float Aux[16][16][32768]);
void Constructor2(float Aux[16][16][32768],float Final[512][32768]);

int main(void)
{


  glp_prob *lp;

  FILE *fp;

  int num_clifford_gates=16384;
  int num_elements_in_matrix = 512;
  int num_rows_in_clifford_group = 524288;

    int i,j,k,l;
    double z;

    float x1,x2,x3,x4,x5,x6,x7,x8,
          x9,x10,x11,x12,x13,x14,x15,x16;


    float (*Final)[32768]   = malloc(512 * 32768 * sizeof *Final);
    float (*Aux)[16][32768]  = malloc(16  * 16 * 32768 * sizeof *Aux);
    float (*Clifford)[16]    = malloc(524288 * 16 * sizeof *Clifford);

    int    *row    = malloc((2*512*num_clifford_gates+1) * sizeof *row);
    int    *column = malloc((2*512*num_clifford_gates+1) * sizeof *column);
    double *values = malloc((2*512*num_clifford_gates+1) * sizeof *values);

  /*-------------------------------------------------------------------------------------------------------------------------*/
  /*Load Clifford matrices. */
    i=0;j=0;
    fp = fopen("Diagonal_Clifford_group_4qubit.txt","r");
    while(fscanf(fp ,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", &x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16 )!=EOF)
      {
        Clifford[i][0]=x1;
        Clifford[i][1]=x2;
        Clifford[i][2]=x3;
        Clifford[i][3]=x4;
        Clifford[i][4]=x5;
        Clifford[i][5]=x6;
        Clifford[i][6]=x7;
        Clifford[i][7]=x8;
        Clifford[i][8]=x9;
        Clifford[i][9]=x10;
        Clifford[i][10]=x11;
        Clifford[i][11]=x12;
        Clifford[i][12]=x13;
        Clifford[i][13]=x14;
        Clifford[i][14]=x15;
        Clifford[i][15]=x16;
        i=i+1;
      }

    fclose(fp);
    i=0;j=0;

 /*Process the Clifford Matrox into the desired form-*----------------------------------------------------------------------------------------*/
    Constructor1(Clifford, Aux);
    Constructor2(Aux ,Final);
  /*-------------------------------------------------------------------------------------------------------------------------*/


  lp = glp_create_prob();
  glp_set_prob_name(lp, "sample");
  glp_set_obj_dir(lp, GLP_MIN);


  glp_add_rows(lp, num_elements_in_matrix); /*rows, elements of the matrix to be decomposed (flattened)*/

   /*Real --------------------------------*/
        for (int i = 1; i <= 256; i++) {
            if (i == 256) {
                glp_set_row_bnds(lp, i, GLP_FX, -1.0, -1.0); // CCCZ: last diagonal element is -1
            } else if ((i - 1) % 17 == 0) {
                glp_set_row_bnds(lp, i, GLP_FX, 1.0, 1.0); // Diagonal elements are 1
            } else {
                glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0); // Off-diagonal elements are 0
            }
        }

    // Set imaginary part constraints for next 256 rows (257-512)
    for (int i = 257; i <= 512; i++) {
        glp_set_row_bnds(lp, i, GLP_FX, 0.0, 0.0);
    }




  // Number of columns is twice the number of clifford groups
  glp_add_cols(lp, (2*num_clifford_gates));
  for ( i = 0; i < (2*num_clifford_gates); i++)
  {
    glp_set_col_bnds(lp,(i+1),GLP_LO,0.0,0.0);
    glp_set_obj_coef(lp, (i+1), 1.0);
  }
  i=0;


  k=1;
  l=1;
  for (j = 1; j <= 512; j++)
  {
    for ( i = 1; i < 2*num_clifford_gates+1; i++)
    {


      row[l]=j;
      column[l]=i;
      values[l]=Final[j-1][i-1];

      l=l+1;



    }

  }

  i=0;j=0;k=0;l=0;
  clock_t start, end;
  double cpu_time_used;

  glp_load_matrix(lp, 2*num_clifford_gates*512, row, column, values);
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


void Constructor1(float Clifford[524288][16], float Aux[16][16][32768])
{
  int i,j,k,l;
  for (i = 0; i < 32768; i++)
  {
    for ( j = 0; j < 16; j++)
    {
      for ( k = 0; k < 16; k++)
      {
        Aux[j][k][i]=Clifford[i*16+j][k];
      }
    }
  }
}


void Constructor2(float Aux[16][16][32768],float Final[512][32768])
{
  int i,j,k,l,count;
  count=0;
  k=0;l=0;
  int num_qubits = 4;
  int matrix_dim = pow(2, num_qubits);
  for ( j = 0; j < 256; j++)/*real part*/
  {
    k = j / matrix_dim;
    l = j % matrix_dim;


    for ( i = 0; i < 16384; i++)
    {
      Final[j][count]=Aux[k][l][2*i];
      Final[j][1+count]=-Aux[k][l][2*i];
      count=count+2;
    }
    count=0;
    k=0;
    l=0;
  }
    for ( j = 256; j < 512; j++)/*Imaginary part*/
  {
    int matrix_index = j - matrix_dim * matrix_dim;
    k = matrix_index / matrix_dim;
    l = matrix_index % matrix_dim;

    count=0;
    for ( i = 0; i < 16384; i++)
    {
      Final[j][count]=Aux[k][l][2*i+1];
      Final[j][1+count]=-Aux[k][l][2*i+1];
      count=count+2;
    }
    count=0;
  }
}
