/*************************************************
 * Laplace Serial C Version
 *
 * Temperature is initially 0.0
 * Boundaries are as follows:
 *
 *      0         T         0
 *   0  +-------------------+  0
 *      |                   |
 *      |                   |
 *      |                   |
 *   T  |                   |  T
 *      |                   |
 *      |                   |
 *      |                   |
 *   0  +-------------------+ 100
 *      0         T        100
 *
 *  John Urbanic, PSC 2014
 *
 ************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

// size of plate
#define COLUMNS    1000
#define ROWS       1000
#define SEC_ROWS   ROWS/4

// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01

double Temperature[ROWS+2][COLUMNS+2];      // temperature grid
double Temperature_last[ROWS+2][COLUMNS+2]; // temperature grid from last iteration

//   helper routines
void initialize();
void track_progress(int iter);


int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    double dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers
    /*
    printf("Maximum iterations [100-4000]?\n");
    scanf("%d", &max_iterations);*/


    // the usual MPI startup routines
    int my_PE_num, numprocs;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);

    // verify only the number of PEs allowed are being used
    if (numprocs != 4 && my_PE_num == 0) {
	printf("Error: incorrect number of processes used! Please only use 4 processes for this problem.");
	return 1;}

    //PE 0 asks for input from user, being the master node, the captain, the teamlead
    if (my_PE_num == 0) {
    	printf("Maximum iterations [100-4000]?\n");
    	scanf("%d", &max_iterations);
    }

    // Now PE 0 has taken the max_iterations from user, broadcast it to other PEs just once.
    // No if clauses wrapped around MPI_Bcast, as MPI_Bcast should be used in every PE, most on receiving end, one on broadcasting end.
    MPI_Bcast(&max_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (my_PE_num == 0) gettimeofday(&start_time,NULL); // Unix timer

    // This initialize() function needs to be modified to take two inputs: numprocs, my_PE_num to assign correct B.C.s for each
    initialize(numprocs, my_PE_num);                   // initialize Temp_last including boundary conditions

    // do until error is minimal (thermal steady-state reached) or until max steps
    while ( dt > MAX_TEMP_ERROR && iteration <= max_iterations ) {

	//printf("In while loop of PE# %d on iteration# %d\n", my_PE_num, iteration);
        // main calculation: average my four neighbors, including the padding boundaries up and down (except for top and bottom pieces)
        if (my_PE_num == 0) {			//top piece
		for(i = 1; i <= SEC_ROWS+1; i++) {
			//if(i % 50 == 1) printf("Entering main calculation for-loop of PE# %d with i=%d at iteration %d\n", my_PE_num, i, iteration);
            		for(j = 1; j <= COLUMNS; j++) {
                		Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            	Temperature_last[i][j+1] + Temperature_last[i][j-1]);
            		}
			//if (j % 50 == 0) printf("Finished calculating Temperature[%d][j] for PE# %d\n", i, my_PE_num);
        	}
	} else if (my_PE_num == numprocs-1) {	//bottom piece
                for(i = 0+SEC_ROWS*my_PE_num; i <= SEC_ROWS+SEC_ROWS*my_PE_num; i++) {
			//if(i % 50 == 0) printf("Entering main calculation for-loop of PE# %d with i=%d at iteration %d\n", my_PE_num, i, iteration);
                        for(j = 1; j <= COLUMNS; j++) {
				//if(j % 50 == 0) printf("About to calculate Temperature[%d][%d] of PE# %d on iteration %d\n", i,j,my_PE_num,iteration);
                        	Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                                Temperature_last[i][j+1] + Temperature_last[i][j-1]);
				//if(j % 50 == 0) printf("Done calculating Temperature[%d][%d] of PE# %d on iteration %d\n", i,j,my_PE_num,iteration);
                        }
			//if (j % 50 == 0) printf("Finished calculating Temperature[%d][j] for PE# %d\n", i, my_PE_num);
                }
	} else {				// other pieces
                for(i = 0+SEC_ROWS*my_PE_num; i <= SEC_ROWS+1+SEC_ROWS*my_PE_num; i++) {
			//if(i % 50 == 0) printf("Entering main calculation for-loop of PE# %d with i=%d at iteration %d\n", my_PE_num, i, iteration);
                        for(j = 1; j <= COLUMNS; j++) {
				//if(j % 50 == 0) printf("About to calculate Temperature[%d][%d] of PE# %d on iteration %d\n", i,j,my_PE_num,iteration);
                        	Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                                Temperature_last[i][j+1] + Temperature_last[i][j-1]);
				//if(j % 50 == 0) printf("Done calculating Temperature[%d][%d] of PE# %d on iteration %d\n", i,j,my_PE_num,iteration);
                        }
			//if (j % 50 == 0) printf("Finished calculating Temperature[%d][j] for PE# %d\n", i, my_PE_num);
                }
	}
	//printf("In while loop of PE# %d, finished main calculation on iteration# %d\n", my_PE_num, iteration);
	// COMMUNICATION PHASE: send and receive padding rows for next iteration
	if (my_PE_num == 0) {
		//printf("Onto Communication task of PE# %d on iteration# %d\n", my_PE_num, iteration);
		MPI_Send(&Temperature[SEC_ROWS+1][1], COLUMNS, MPI_DOUBLE, 1, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, 1, iteration);
		MPI_Recv(&Temperature[SEC_ROWS+1][1], COLUMNS, MPI_DOUBLE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", 1, my_PE_num, iteration);
	} else if (my_PE_num == numprocs-1) {
		//printf("Onto Communication task of PE# %d on iteration# %d\n", my_PE_num, iteration);
		MPI_Send(&Temperature[SEC_ROWS*my_PE_num][1], COLUMNS, MPI_DOUBLE, numprocs-2, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, numprocs-2, iteration);
		MPI_Recv(&Temperature[SEC_ROWS*my_PE_num][1], COLUMNS, MPI_DOUBLE, numprocs-2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", numprocs-2, my_PE_num, iteration);
	} else {
		//printf("Onto Communication task of PE# %d on iteration# %d\n", my_PE_num, iteration);
		MPI_Send(&Temperature[SEC_ROWS*my_PE_num][1], COLUMNS, MPI_DOUBLE, my_PE_num-1, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, my_PE_num-1, iteration);
		MPI_Recv(&Temperature[SEC_ROWS*my_PE_num][1], COLUMNS, MPI_DOUBLE, my_PE_num-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num-1, my_PE_num, iteration);
		MPI_Send(&Temperature[SEC_ROWS*my_PE_num+SEC_ROWS+1][1], COLUMNS, MPI_DOUBLE, my_PE_num+1, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, my_PE_num+1, iteration);
		MPI_Recv(&Temperature[SEC_ROWS*my_PE_num+SEC_ROWS+1][1], COLUMNS, MPI_DOUBLE, my_PE_num+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num+1, my_PE_num, iteration);
	}
	//printf("In while loop of PE# %d, finished communication phase on iteration# %d\n", my_PE_num, iteration);
	// set barrier here to make sure all current time-step calculations are synced before moving to next iteration
	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Went past MPI Barrier on PE# %d\n", my_PE_num);
        dt = 0.0; // reset largest temperature change

        // copy grid to old grid for next iteration, and find latest dt to see if steady-state is reached
        for(i = SEC_ROWS*my_PE_num+1; i <= SEC_ROWS*my_PE_num+SEC_ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
	      dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
	      Temperature_last[i][j] = Temperature[i][j];
//	      if (i == 995 && j == 995) printf("Temperature_last[995][995]=%f and Temperature_last[1001][1001] = %f\n",Temperature_last[995][995],Temperature_last[1001][1001]);
            }
        }
	//printf("Finished updating Temperature_last grid for next iteration on PE# %d\n", my_PE_num);

        // periodically print test values
        if((iteration % 10) == 0 && my_PE_num == 0) {
 	    track_progress(iteration);
        }

	iteration++;
    }

    if (my_PE_num == 0) {
    	gettimeofday(&stop_time,NULL);
    	timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    	printf("\nMax error at iteration %d was %f\n", iteration-1, dt);
    	printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    }

}


// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize(int numprocs, int my_PE_num){

    int i,j;
    printf("Now inside initialize function of PE# %d\n", my_PE_num);
    // step 1: initialize all to 0, and set left side and right side
    for(i = SEC_ROWS*my_PE_num; i <= SEC_ROWS*my_PE_num+SEC_ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        }
	Temperature_last[i][COLUMNS+1] = my_PE_num/numprocs*100.0 + (100.0/ROWS)*i;
	if(i % 50 == 0) printf("Temperature_last[%d][%d] initialized as %f\n",i,COLUMNS+1,Temperature_last[i][COLUMNS+1]);
    }

    printf("Finished step 1, now on the step 2 of initialize function of PE# %d\n", my_PE_num);
    // step 2: set top side and bottom side
    if (my_PE_num == numprocs-1){
    	for(j = 0; j <= COLUMNS+1; j++) {
        	Temperature_last[SEC_ROWS*my_PE_num][j] = 0.0;
        	Temperature_last[SEC_ROWS*my_PE_num+SEC_ROWS+1][j] = (100.0/COLUMNS)*j;
		if(j % 50 == 0) printf("Temperature_last[%d][%d] initialized as %f\n",COLUMNS+1,j,Temperature_last[COLUMNS+1][j]);
    	}
    } else {
	for(j = 0; j<= COLUMNS+1; j++) {
		Temperature_last[SEC_ROWS*my_PE_num][j] = 0.0;
		Temperature_last[SEC_ROWS*my_PE_num+SEC_ROWS+1][j] = 0.0;
	}
    }
    printf("Finished step 2 of initialization function of PE# %d\n", my_PE_num);
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {

    int i;

    printf("---------- Iteration number: %d ------------\n", iteration);
    for(i = ROWS-5; i <= ROWS; i++) {
        printf("[%d,%d]: %5.2f  ", i, i, Temperature_last[i][i]);
    }
    printf("\n");
}
