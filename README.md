# HPC_Exercise
HPC Expercises from XSEDE HPC workshop

## MPI Terminal Commands
To compile a c-code (hello.c for example) into executable:
```shell script
mpicc hello.c
```
To run the compiled code (a.out for example) using MPI:
```shell script
mpirun -n 8 a.out
```
An example of send and receive is made in hello_sendrecv.c, to run it with two processes:
```shell script
mpicc hello_sendrecv.c
mpirun -n 2 a.out
```
