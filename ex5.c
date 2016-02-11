#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int max(int a, int b){
  if(a<b) return b;
  else return a;
}

int main (int argc, char** argv){
    int np, rank;
    int rcs = 16; //num rows/cols of A.To make it easier: rcs divisible by np
    int row, column;
    double *A = NULL;
    A = malloc(rcs*rcs*sizeof(double));
    FILE *fp;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //read matrix A
    if (rank==0){
      fp = fopen("matrika_eva.out", "r");
      for (row = 0; row < rcs; row++){
        for (column = 0; column < rcs; column++)
          fscanf(fp, "%lf", &A[row*rcs+column]);
      }
    }
    //scatter rows to all. count = rows per rank
    int count = rcs/np; //dont forget! no remainder allowed!
    double *recv_A = malloc(count*rcs*sizeof(double));

    if(rank==0) printf("count=%d\n", count);

    MPI_Scatter(A, count*rcs, MPI_DOUBLE, recv_A,
                count*rcs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Scatter done.\n");

    //start algorithm
    double temp_l; //, akk;
    double * temp_Akj = malloc(rcs*sizeof(double));
    int root;
    for(int k=0; k<rcs-1; k++){
        root = k/count;
        //akk = recv_A[(k%count)*(rcs+1)];
        printf("rank=%d,   k=%d\n", rank, k);
        if(rank==root){
//          printf("k=%d, rank: %d, writing Akj...\n", k, rank);
          for(int ind=0; ind<rcs; ind++){
            temp_Akj[ind] = recv_A[(k%count)*rcs + ind];
          }
        }
        MPI_Bcast(&temp_Akj, rcs, MPI_DOUBLE, root, MPI_COMM_WORLD);
        printf("k=%d, broadcast done, rank %d \n", k, rank);
        //MPI_Bcast(&akk, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        for(int i=max(k+1, rank*count); i<(rank+1)*count; i++){
          temp_l = recv_A[(i%count)*rcs + k]/temp_Akj[k]; //akk;
          for(int j=k; j<rcs; j++){
              recv_A[(i%count)*rcs + j] -= temp_l*temp_Akj[j];

          }
        }
    }

    //Gather the results
    MPI_Gather(recv_A, count*rcs, MPI_DOUBLE, A,
                count*rcs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //print result
    if(rank==0){
      for (row = 0; row < rcs; row++){
        for (column = 0; column < rcs; column++){
          printf("%lf ", A[row*rcs+column]);
        }
        printf("\n");
      }
    }

    MPI_Finalize();

    return 0;
}
