#include <stdio.h>
#include "mpi.h"

main(int argc, char** argv){

  int n, my_PE_num, numbertosend, numbertoreceive;
  MPI_Status status;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);

  numbertosend = my_PE_num * 2;
  n = 8;
/*
  if(my_PE_num == 0) {
	MPI_Recv( &numbertoreceive, 1, MPI_INT, 7, 666, MPI_COMM_WORLD, &status);
  } else MPI_Recv( &numbertoreceive, 1, MPI_INT, my_PE_num-1, 666, MPI_COMM_WORLD, &status);

  MPI_Bcast( &n, 1, MPI_INT, my_PE_num, 666, MPI_COMM_WORLD);
*/

  if (my_PE_num == 7)
	MPI_Send( &numbertosend, 1, MPI_INT, 0, 666, MPI_COMM_WORLD);
  else
	MPI_Send( &numbertosend, 1, MPI_INT, my_PE_num+1, 666, MPI_COMM_WORLD);

  MPI_Recv( &numbertoreceive, 1, MPI_INT, MPI_ANY_SOURCE, 666, MPI_COMM_WORLD, &status);

  printf("Hello from %d and I received %d.\n", my_PE_num, numbertoreceive);

  MPI_Finalize();

}
