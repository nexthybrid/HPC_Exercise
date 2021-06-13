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

// !!!IMPORTANT: EACH WORKER STORES 252x1002 GRID DATA FROM THE BEGINNING, AND COMMUNICATE EACH ITERATION
double Temperature[SEC_ROWS+2][COLUMNS+2];      // temperature grid 252x1002, a temporary variable
double Temperature_last[SEC_ROWS+2][COLUMNS+2]; // temperature grid from last iteration 252x1002

//   helper routines declaration
void initialize(int, int);
void track_progress(int iteration);
void output(int my_PE_num, int iteration);
void output_full(int my_PE_num, int iteration);

int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    double dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

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
    MPI_Barrier(MPI_COMM_WORLD);			// use barrier to make sure all grid is initialized before taking the output(0)
    output_full(my_PE_num, 0);

    // do until error is minimal (thermal steady-state reached) or until max steps
    while ( dt > MAX_TEMP_ERROR && iteration <= max_iterations ) {

	//printf("In while loop of PE# %d on iteration# %d\n", my_PE_num, iteration);
        // main calculation: average my four neighbors, including the padding boundaries up and down (except for top and bottom pieces)

	for(i = 1; i <= SEC_ROWS; i++) {
		//if(i % 50 == 1) printf("Entering main calculation for-loop of PE# %d with i=%d at iteration %d\n", my_PE_num, i, iteration);
            	for(j = 1; j <= COLUMNS; j++) {
                	Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            	    Temperature_last[i][j+1] + Temperature_last[i][j-1]);
		}
	}
	//printf("In while loop of PE# %d, finished main calculation on iteration# %d\n", my_PE_num, iteration);

	// COMMUNICATION PHASE: send and receive padding rows for next iteration
	// Sending and receiving order (after main calculation of Temperature, which is a temporary variable):
	//
	// 	PE0's row 250 	in Temperature 	-> 	PE1's row 0 	in Temperature_last
	//
	//	PE1's row 1 	in Temperature	-> 	PE0's row 251 	in Temperature_last
	//	PE1's row 250	in Temperature	->	PE2's row 0	in Temperature_last
	//
	//	PE2's row 1	in Temperature	->	PE1's row 251	in Temperature_last
	//	PE2's row 250	in Temperature	->	PE3's row 0	in Temperature_last
	//
	//	PE3's row 1	in Temperature	->	PE2's row 251	in Temperature_last
	//
	//	(!NOTE: all values to send are about Temperature array, not Temperature_last array)
	//	(!NOTE: all values to receive are about Temperature_last array, not Temperature array)
	if (my_PE_num == 0) {
		//printf("Onto Communication task of PE# %d on iteration# %d\n", my_PE_num, iteration);
		MPI_Send(&Temperature[SEC_ROWS][1], COLUMNS, MPI_DOUBLE, 1, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, 1, iteration);
		MPI_Recv(&Temperature_last[SEC_ROWS+1][1], COLUMNS, MPI_DOUBLE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", 1, my_PE_num, iteration);
	} else if (my_PE_num == numprocs-1) {
		//printf("Onto Communication task of PE# %d on iteration# %d\n", my_PE_num, iteration);
		MPI_Send(&Temperature[1][1], COLUMNS, MPI_DOUBLE, numprocs-2, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, numprocs-2, iteration);
		MPI_Recv(&Temperature_last[0][1], COLUMNS, MPI_DOUBLE, numprocs-2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", numprocs-2, my_PE_num, iteration);
	} else {
		//printf("Onto Communication task of PE# %d on iteration# %d\n", my_PE_num, iteration);
		MPI_Send(&Temperature[1][1], COLUMNS, MPI_DOUBLE, my_PE_num-1, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, my_PE_num-1, iteration);
		MPI_Recv(&Temperature_last[0][1], COLUMNS, MPI_DOUBLE, my_PE_num-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num-1, my_PE_num, iteration);
		MPI_Send(&Temperature[SEC_ROWS][1], COLUMNS, MPI_DOUBLE, my_PE_num+1, 666, MPI_COMM_WORLD);
		//printf("Send from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num, my_PE_num+1, iteration);
		MPI_Recv(&Temperature_last[SEC_ROWS+1][1], COLUMNS, MPI_DOUBLE, my_PE_num+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//printf("Receive from PE# %d to PE# %d on iteration# %d done!\n", my_PE_num+1, my_PE_num, iteration);
	}
	//printf("In while loop of PE# %d, finished communication phase on iteration# %d\n", my_PE_num, iteration);
	// set barrier here to make sure all current time-step calculations are synced before moving to next iteration
	MPI_Barrier(MPI_COMM_WORLD);

	//printf("Went past MPI Barrier on PE# %d\n", my_PE_num);
        dt = 0.0; // reset largest temperature change

        // copy grid to old grid for next iteration, and find latest dt to see if steady-state is reached
        for(i = 1; i <= SEC_ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
	      dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
	      Temperature_last[i][j] = Temperature[i][j];
	      // if (i == 995 && j == 995) printf("Temperature_last[995][995]=%f and Temperature_last[1001][1001] = %f\n",Temperature_last[995][995],Temperature_last[1001][1001]);
            }
        }
	//printf("Finished updating Temperature_last grid for next iteration on PE# %d\n", my_PE_num);

        // periodically print test values
        if((iteration % 500) == 0 && my_PE_num == 3) {
 	    	track_progress(iteration);
		fflush(stdout);
        }

	// periodically print temperature grid
	if(iteration % 500 == 0) {
		output_full(my_PE_num, iteration);
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


// initialize plate and boundary conditions on the 252x1002 grid for each PE
// Temp_last is used to to start first iteration
void initialize(int numprocs, int my_PE_num){

    int i,j;
    printf("Now inside initialize function of PE# %d\n", my_PE_num);
    // step 1: initialize all 252x1002 to 0 for cleanup, and set left side and right side
    for(i = 0; i<= SEC_ROWS+1; i++) {	// now every PE gets to initialize the whole 252x1002
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        };
	Temperature_last[i][COLUMNS+1] = my_PE_num*100.0/numprocs + (100.0/ROWS)*i; // IMPORTANT: the order of operation matters!
	// In the above line, if my_PE_num/numprocs*100.0 is used, the result will always be zero due to int/int produces an int with rounding!
	if(i % 50 == 0) printf("Temperature_last[%d][%d] of PE# %d initialized as %f\n",i,COLUMNS+1, my_PE_num, Temperature_last[i][COLUMNS+1]);
    }
    //printf("Finished step 1, now on the step 2 of initialize function of PE# %d\n", my_PE_num);

    // step 2: set top side and bottom side, WARNING: DON'T RUIN THE LEFT-RIGHT INITILIZATION ABOVE
    if (my_PE_num == numprocs-1){	// only the bottom piece has a heating element on horizontal boundary
    	for(j = 0; j <= COLUMNS+1; j++) {
        	Temperature_last[SEC_ROWS+1][j] = (100.0/COLUMNS)*j;
		if(j % 200 == 0) printf("Temperature_last[%d][%d] initialized as %f\n",SEC_ROWS+1,j,Temperature_last[SEC_ROWS+1][j]);
    	}
    }
    printf("Initialization Done on PE# %d\n", my_PE_num);
}

// print diagonal in bottom right corner of the 252x1002 piece
void track_progress(int iteration) {
    int i;
    printf("---------- Iteration number: %d ------------\n", iteration);
    for(i = 0; i <= 5; i++) {
        printf("[%d,%d]: %5.2f  ", SEC_ROWS-i, COLUMNS-i, Temperature_last[SEC_ROWS-i][COLUMNS-i]);
    };
    printf("\n");
    printf("Finished track_progress(%d)", iteration);
}

// output the semi-full 1000x1002 grid of temperatures. Using MPI_Barrier to make sure all PEs are caught up to print
void output(int my_PE_num, int iteration) {
    FILE* fp;
    char filename[50];
    sprintf(filename,"output%d.csv",iteration);

	for (int pe = 0; pe <= 3; pe++) {
		if(my_PE_num == pe) {
			fp = fopen(filename, "a");	// "a" for append to a file, "w" for write / rewrite a file

    			for(int i = 1; i <= SEC_ROWS; i++){	// for each PE, print 250x1002 grid, includ. left and right padding
				for(int j = 0; j <= COLUMNS; j++){
					fprintf(fp, "%5.2f,",Temperature_last[i][j]);
				};
				fprintf(fp, "%5.2f\n",Temperature_last[i][COLUMNS+1]); // use \n to start a new row in csv file
    			};
    			fflush(fp);
    			fclose(fp);
		};
		MPI_Barrier(MPI_COMM_WORLD);	// use barrier inside the pe for loop to stop all PEs until all PEs have printed their share
	}
}

// output the full 1002x1002 grid of temperatures including the horizontal and vertical padding boundaries.
void output_full(int my_PE_num, int iteration) {
	FILE* fp2;
	char filename[50];
	sprintf(filename, "outputf%d.csv", iteration);

        for (int pe = 0; pe <= 3; pe++) {
                if(my_PE_num == pe) {
                        fp2 = fopen(filename, "a");      // "a" for append to a file, "w" for write / rewrite a file

			//for PE = 0 add the top row with 1x1002 elements
			if(my_PE_num == 0) {
				for(int jj = 0; jj <= COLUMNS; jj++) {
					fprintf(fp2, "%5.2f,", Temperature_last[0][jj]);
				};
				fprintf(fp2, "5.2f\n", Temperature_last[0][COLUMNS+1]);
			};

			// for any PE print the 250x1002 grid elements
                        for(int i = 1; i <= SEC_ROWS; i++){     // for each PE, print 250x1002 grid, includ. left and right padding
                                for(int j = 0; j <= COLUMNS; j++){
                                        fprintf(fp2, "%5.2f,",Temperature_last[i][j]);
                                };
                                fprintf(fp2, "%5.2f\n",Temperature_last[i][COLUMNS+1]); // use \n to start a new row in csv file
                        };

			// for PE = 3 add the bottom padding row with 1x1002 elements
			if(my_PE_num == 3) {
				for(int jjj = 0; jjj <= COLUMNS; jjj++) {
					fprintf(fp2, "%5.2f,", Temperature_last[SEC_ROWS+1][jjj]);
					//printf("%5.2f",Temperature_last[SEC_ROWS+1][jjj]);	// for debug
				};
				fprintf(fp2, "5.2f\n", Temperature_last[SEC_ROWS+1][COLUMNS+1]);
			};
                        fflush(fp2);
                        fclose(fp2);
                };
                MPI_Barrier(MPI_COMM_WORLD);    // use barrier inside the pe for loop,
		// such that it blocks progress at each loop: PE = 0, PE = 1, PE = 2, PE = 3,
		// to stop all PEs until all PEs have printed their share in the correct order.
        }
}
