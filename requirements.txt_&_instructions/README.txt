The Project is a parallelization exercise for Nonlinear Optimization using the algorithm of Hooke and Jeeves 12 February 1994.
Writen 100% in C using mainly OpenMp and MPI interfaces.
	--Patra_Ceid::Parallel_Programming_Project_2022--

Instructions:

	4 parallelization methods are used : OpenMP, OpenMP Tasks, MPI, MPI+OpenMP(hybrid model)
	Along side those the sequential code is accounted for with also 2 bash scripts, a Makefile 
	and a python script.
	
	Creation of a python Virtual environment is advised with the installed requirements found
	in this folder. Use of Linux and bash shell is also highly advised or small changes to the 
	scipts must be made.
	
	multipleExecutionScript.sh is used with 3 command line arguments - parameters:
		1) # of trials -> problem dimension
		2) # of max threads the executables will use
		3) # of optimizations used with gcc
	The generated txts are the Results of all the corresponding parallelization methods,
	till the number of threads specified in argument 2. At last the python script 
	visualizes the results with a Strong scaling Speedup(p) and efficiency.
