# HPC_Exercise
HPC Expercises from XSEDE HPC workshop

## MPI Terminal Commands
To compile a c-code (hello.c for example) into executable:
```shell script
mpicc hello.c
```
To run the compiled code (a.out for example) using MPI with 8 processes:
```shell script
mpirun -n 8 a.out
```
An example of send and receive is made in hello_sendrecv.c, to run it with two processes:
```shell script
mpicc hello_sendrecv.c
mpirun -n 2 a.out
```

## File Transfer Commands
To transfer data from a bridges2 HPC node to personal PC, on personal PC terminal type:
```shell script
scp zhao1991@bridges2.psc.edu:~/HPC_Exercise/Exercises/MPI/outputf0.csv c:\github\HPC_Exercise\Exercises\MPI\
```
The syntax format is `scp [xsede_user_name]@bridges2.psc.edu:[bridges2 file path] [personal pc file path]`. You'll be asked to enter PSC password to complete the transfer.

## Debugging Tips for MPI Programming
1. Program keeps running, while print statements got stuck. 
Possible reason: 1) basic c grammar error, like 
```C
if(int i = 0; i < 5; i++)	// should be ==
```
2) running function/code intended for all PEs(usually MPI_Barrier is placed inside the function/code) in a local block for one PE. For example, if some code is intended for all PEs to go through, but you run in instead only for PE0:
```C
void func_for_all_pe(int pe_num) {
	// do something in a for loop in pe = 0,1,2... order, each pe may have slightly different behavior
	for(int pe = 0; pe < numprocs; pe++) {
		if(pe_num == pe) {
			// PE-specific behavior, i.e. writing this PE's calculated data to file
		};
		MPI_Barrier(MPI_COMM_WORLD); // Blocking program across all workers untill all workers get here
	}
}

int main(int argc, char *argv[]) {
	// call func_for_all_pe() in a local bock just for PE = 0
	if(pe_num == 0) {
		func_for_all_pe(pe_num);
	}
	// rest of the code
}
```
In this case, you are running `func_for_all_pe()` only for PE0, so all other PEs don't get to access `func_for_all_pe`, so the `MPI_Barrier(MPI_COMM_WORLD)` will never fulfuill and let pass, causing the program to stall.