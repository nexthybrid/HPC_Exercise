% visualize Laplace heat transfer plat

% The scp command for porting data out of bridge2 node
% scp zhao1991@bridges2.psc.edu:~/HPC_Exercise/Exercises/MPI/outputf0.csv c:\github\HPC_Exercise\Exercises\MPI\

tempData = readmatrix("outputf500.csv");

surf(tempData,'EdgeColor','none')
xlabel("x")
ylabel("y")
zlabel("temperature")