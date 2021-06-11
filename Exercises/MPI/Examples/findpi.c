#include <mpi.h>
#include <math.h>
#include <stdio.h>

int main( int argc, char **argv ){
 	int n, my_pe_num, numprocs, index;
	float mypi, pi, h, x, start, end;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); //this automatically decides the number of processes from user input
	MPI_Comm_rank(MPI_COMM_WORLD, &my_pe_num);	//assign the rank # to &my_pe_num

	if ( my_pe_num == 0 ){
		printf("How many intervals? ");
		scanf("%d", &n);
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);	//both master and worker will use the same MPI_Bcast() command, but only master(process 0) will be broadcasting.
	// MPI_Bcast is also an (implicit) barrier.

	mypi = 0;
	h = (float) 2/n;	//size of each slice
	start = (my_pe_num*2/numprocs)-1;	//slices for this PE
	end = ((my_pe_num+1)*2/numprocs)-1;

	for (x = start; x < end; x = x+h)
		mypi = mypi + h * 2 * sqrt(1-x*x);

	MPI_Reduce(&mypi, &pi, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	if ( my_pe_num == 0 ){
		printf("Pi is approximately %f\n", pi);
		printf("Error is %f\n", pi-3.141592653589779323846);
	}

	MPI_Finalize();
}
