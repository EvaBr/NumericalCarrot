#include <stdio.h>
#include <stdlib.h>

int main (int argc, char** argv){
    int rcs = 16;
    int row, column;
    double temp_l;
    double *A = NULL;
    A = malloc(rcs*rcs*sizeof(double));
    FILE *fp;
    fp = fopen("matrika_eva.out", "r");

    //read matrix A
    for (row = 0; row < rcs; row++){
      for (column = 0; column < rcs; column++){
        fscanf(fp, "%lf,", &A[row*rcs+column]);
        printf("%lf ", A[row*rcs+column]);
      }
      printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    //start algorithm
    for(int k=0; k<rcs-1; k++){
        for(int i=k+1; i<rcs; i++){
          temp_l = A[i*rcs + k]/A[k*rcs+k];
          for(int j=k; j<rcs; j++){
              A[i*rcs + j] -= temp_l*A[k*rcs+j];
          }
        }
    }

    //print result
    for (row = 0; row < rcs; row++){
      for (column = 0; column < rcs; column++){
        printf("%lf ", A[row*rcs+column]);
      }
      printf("\n");
    }
}
