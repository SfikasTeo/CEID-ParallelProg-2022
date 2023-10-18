/* Nonlinear Optimization using the algorithm of Hooke and Jeeves  */
/*	12 February 1994	author: Mark G. Johnson 	   */


/* Find a point X where the nonlinear function f(X) has a local    */
/* minimum.  X is an n-vector and f(X) is a scalar.  In mathe-	   */
/* matical notation  f: R^n -> R^1.  The objective function f()    */
/* is not required to be continuous.  Nor does f() need to be	   */
/* differentiable.  The program does not use or require 	   */
/* derivatives of f().						   					*/

/* The software user supplies three things: a subroutine that	   */
/* computes f(X), an initial "starting guess" of the minimum point */
/* X, and values for the algorithm convergence parameters.  Then   */
/* the program searches for a local minimum, beginning from the    */
/* starting guess, using the Direct Search algorithm of Hooke and  */
/* Jeeves.							   */

/* This C program is adapted from the Algol pseudocode found in    */
/* "Algorithm 178: Direct Search" by Arthur F. Kaupe Jr., Commun-  */
/* ications of the ACM, Vol 6. p.313 (June 1963).  It includes the */
/* improvements suggested by Bell and Pike (CACM v.9, p. 684, Sept */
/* 1966) and those of Tomlin and Smith, "Remark on Algorithm 178"  */
/* (CACM v.12).  The original paper, which I don't recommend as    */
/* highly as the one by A. Kaupe, is:  R. Hooke and T. A. Jeeves,  */
/* "Direct Search Solution of Numerical and Statistical Problems", */
/* Journal of the ACM, Vol. 8, April 1961, pp. 212-229. 	   */

/* Calling sequence:						   */
/*  int hooke(nvars, startpt, endpt, rho, epsilon, itermax)	   */
/*								   */
/*     nvars	   {an integer}  This is the number of dimensions  */
/*		   in the domain of f().  It is the number of	   */
/*		   coordinates of the starting point (and the	   */
/*		   minimum point.)				   */
/*     startpt	   {an array of doubles}  This is the user-	   */
/*		   supplied guess at the minimum.		   */
/*     endpt	   {an array of doubles}  This is the location of  */
/*		   the local minimum, calculated by the program    */
/*     rho	   {a double}  This is a user-supplied convergence */
/*		   parameter (more detail below), which should be  */
/*		   set to a value between 0.0 and 1.0.	Larger	   */
/*		   values of rho give greater probability of	   */
/*		   convergence on highly nonlinear functions, at a */
/*		   cost of more function evaluations.  Smaller	   */
/*		   values of rho reduces the number of evaluations */
/*		   (and the program running time), but increases   */
/*		   the risk of nonconvergence.	See below.	   */
/*     epsilon	   {a double}  This is the criterion for halting   */
/*		   the search for a minimum.  When the algorithm   */
/*		   begins to make less and less progress on each   */
/*		   iteration, it checks the halting criterion: if  */
/*		   the stepsize is below epsilon, terminate the    */
/*		   iteration and return the current best estimate  */
/*		   of the minimum.  Larger values of epsilon (such */
/*		   as 1.0e-4) give quicker running time, but a	   */
/*		   less accurate estimate of the minimum.  Smaller */
/*		   values of epsilon (such as 1.0e-7) give longer  */
/*		   running time, but a more accurate estimate of   */
/*		   the minimum. 				   */
/*     itermax	   {an integer}  A second, rarely used, halting    */
/*		   criterion.  If the algorithm uses >= itermax    */
/*		   iterations, halt.				   */


/* The user-supplied objective function f(x,n) should return a C   */
/* "double".  Its  arguments are  x -- an array of doubles, and    */
/* n -- an integer.  x is the point at which f(x) should be	   */
/* evaluated, and n is the number of coordinates of x.	That is,   */
/* n is the number of coefficients being fitted.		   */

/* rho, the algorithm convergence control			   */
/*	The algorithm works by taking "steps" from one estimate of */
/*    a minimum, to another (hopefully better) estimate.  Taking   */
/*    big steps gets to the minimum more quickly, at the risk of   */
/*    "stepping right over" an excellent point.  The stepsize is   */
/*    controlled by a user supplied parameter called rho.  At each */
/*    iteration, the stepsize is multiplied by rho  (0 < rho < 1), */
/*    so the stepsize is successively reduced.			   */
/*	Small values of rho correspond to big stepsize changes,    */
/*    which make the algorithm run more quickly.  However, there   */
/*    is a chance (especially with highly nonlinear functions)	   */
/*    that these big changes will accidentally overlook a	   */
/*    promising search vector, leading to nonconvergence.	   */
/*	Large values of rho correspond to small stepsize changes,  */
/*    which force the algorithm to carefully examine nearby points */
/*    instead of optimistically forging ahead.	This improves the  */
/*    probability of convergence.				   */
/*	The stepsize is reduced until it is equal to (or smaller   */
/*    than) epsilon.  So the number of iterations performed by	   */
/*    Hooke-Jeeves is determined by rho and epsilon:		   */
/*	    rho**(number_of_iterations) = epsilon		   */
/*	In general it is a good idea to set rho to an aggressively */
/*    small value like 0.5 (hoping for fast convergence).  Then,   */
/*    if the user suspects that the reported minimum is incorrect  */
/*    (or perhaps not accurate enough), the program can be run	   */
/*    again with a larger value of rho such as 0.85, using the	   */
/*    result of the first minimization as the starting guess to    */
/*    begin the second minimization.				   */

/* Normal use: (1) Code your function f() in the C language	   */
/*	       (2) Install your starting guess {or read it in}	   */
/*	       (3) Run the program				   */
/*	       (4) {for the skeptical}: Use the computed minimum   */
/*		      as the starting point for another run	   */

/* Data Fitting:						   */
/*	Code your function f() to be the sum of the squares of the */
/*	errors (differences) between the computed values and the   */
/*	measured values.  Then minimize f() using Hooke-Jeeves.    */
/*	EXAMPLE: you have 20 datapoints (ti, yi) and you want to   */
/*	find A,B,C such that  (A*t*t) + (B*exp(t)) + (C*tan(t))    */
/*	fits the data as closely as possible.  Then f() is just    */
/*	f(x) = SUM (measured_y[i] - ((A*t[i]*t[i]) + (B*exp(t[i])) */
/*				  + (C*tan(t[i]))))^2		   */
/*	where x[] is a 3-vector consisting of {A, B, C}.	   */

/*								   */
/*  The author of this software is M.G. Johnson.		   */
/*  Permission to use, copy, modify, and distribute this software  */
/*  for any purpose without fee is hereby granted, provided that   */
/*  this entire notice is included in all copies of any software   */
/*  which is or includes a copy or modification of this software   */
/*  and in all copies of the supporting documentation for such	   */
/*  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT    */
/*  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE   */
/*  AUTHOR NOR AT&T MAKE ANY REPRESENTATION OR WARRANTY OF ANY	   */
/*  KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS    */
/*  FITNESS FOR ANY PARTICULAR PURPOSE. 			   				*/
/*								   									*/

//mpiexec -n # executable #arguments
#define __STDC_WANT_LIB_EXT2__ 1  //Define you want TR 24731-2:2010 extensions -- for gcc -- as printf

#include <stdio.h> // i/o
#include <string.h>
#include <stdlib.h> //Memory management
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

#define MAXVARS		(32)	/* max # of variables	     */
#define RHO_BEGIN	(0.9)	/* stepsize geometric shrink */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */

/* Rosenbrock classic parabolic valley ("banana") function */
double f(double *x, int n, unsigned long *funevals)
{
    double fv;
    int i;

    (*funevals)++;
    fv = 0.0;
    for (i=0; i<n-1; i++)   /* rosenbrock */
        fv = fv + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);

    return fv;
}

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double delta[MAXVARS], double point[MAXVARS], double prevbest, int nvars, unsigned long *funevals )
{
	double z[MAXVARS];
	double minf, ftmp;
	int i;
	minf = prevbest;
	for (i = 0; i < nvars; i++)
		z[i] = point[i];
	for (i = 0; i < nvars; i++) {
		z[i] = point[i] + delta[i];
		ftmp = f(z, nvars, funevals);
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i] = point[i] + delta[i];
			ftmp = f(z, nvars,funevals);
			if (ftmp < minf)
				minf = ftmp;
			else
				z[i] = point[i];
		}
	}
	for (i = 0; i < nvars; i++)
		point[i] = z[i];

	return (minf);
}

int hooke(int nvars, double startpt[MAXVARS], double endpt[MAXVARS], double rho, double epsilon, int itermax ,unsigned long *funevals )
{
	double delta[MAXVARS];
	double newf, fbefore, steplength, tmp;
	double xbefore[MAXVARS], newx[MAXVARS];
	int i, j, keep;
	int iters, iadj;

	for (i = 0; i < nvars; i++) {
		newx[i] = xbefore[i] = startpt[i];
		delta[i] = fabs(startpt[i] * rho);
		if (delta[i] == 0.0)
			delta[i] = rho;
	}
	
	iadj = 0;
	steplength = rho;
	iters = 0;
	fbefore = f(newx, nvars,funevals);
	newf = fbefore;
	while ((iters < itermax) && (steplength > epsilon)) {
		iters++;
		iadj++;
#if DEBUG
		printf("\nAfter %5d funevals, f(x) =  %.4le at\n", *l_funevals, fbefore);
		for (j = 0; j < nvars; j++)
			printf("   x[%2d] = %.4le\n", j, xbefore[j]);
#endif
		/* find best new point, one coord at a time */
		for (i = 0; i < nvars; i++) {
			newx[i] = xbefore[i];
		}
		newf = best_nearby(delta, newx, fbefore, nvars, funevals);
		/* if we made some improvements, pursue that direction */
		keep = 1;
		while ((newf < fbefore) && (keep == 1)) {
			iadj = 0;
			for (i = 0; i < nvars; i++) {
				/* firstly, arrange the sign of delta[] */
				if (newx[i] <= xbefore[i])
					delta[i] = 0.0 - fabs(delta[i]);
				else
					delta[i] = fabs(delta[i]);
				/* now, move further in this direction */
				tmp = xbefore[i];
				xbefore[i] = newx[i];
				newx[i] = newx[i] + newx[i] - tmp;
			}
			fbefore = newf;
			newf = best_nearby(delta, newx, fbefore, nvars, funevals);
			/* if the further (optimistic) move was bad.... */
			if (newf >= fbefore)
				break;

			/* make sure that the differences between the new */
			/* and the old points are due to actual */
			/* displacements; beware of roundoff errors that */
			/* might cause newf < fbefore */
			keep = 0;
			for (i = 0; i < nvars; i++) {
				keep = 1;
				if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i])))
					break;
				else
					keep = 0;
			}
		}
		if ((steplength >= epsilon) && (newf >= fbefore)) {
			steplength = steplength * rho;
			for (i = 0; i < nvars; i++) {
				delta[i] *= rho;
			}
		}
	}
	for (i = 0; i < nvars; i++)
		endpt[i] = xbefore[i];

	return (iters);
}

//custom reduction fuction
void cust_min(double* invec, double* inoutvec, int *len, MPI_Datatype *dtype){

	inoutvec[0] += invec[0];
	if( invec[1] < inoutvec[1] ) inoutvec[1] = invec[1];
}

int main(int argc, char **argv)
{
	MPI_Init(&argc , &argv); //Initializes the MPI environment

        //Mpi basics
		int rank, size;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank); //initialize the rank using the default communicator
		MPI_Comm_size(MPI_COMM_WORLD,&size); //initialzie the size using the default communicator
	
		__attribute__((aligned(64))) double startpt[MAXVARS];
		__attribute__((aligned(64))) double endpt[MAXVARS];

		int itermax = IMAX;
		double rho = RHO_BEGIN;
		double epsilon = EPSMIN;
		int ntrials = 64*1024;		/* number of trials */
		int nvars = 32;			/* number of variables (problem dimension) */
			
		double fx;
		int i, jj, trial; //i is the global counter 
		double t0, t1;

		//initialization of needed variables
		double best_fx = 1e10;
		double best_pt[MAXVARS];
		
		int best_trial = -1;
		int best_jj = -1;
        	unsigned long funevals = 0;
		long ltimeSeed = time(0);
		unsigned short randBuffer[3] = {0,0,ltimeSeed+rank};

		if (argc > 1) //set the ntrials as the first execution argument
		{
			ntrials = atoi(argv[1]);
		}

		//start the clock
		MPI_Barrier(MPI_COMM_WORLD); //Barrier in order to measure time correctly
		t0 = MPI_Wtime(); 
				
		for (trial = rank; trial < ntrials; trial+=size) {

			//Initialization of startpt randomly with erand
			for (i = 0; i < nvars; i++) {
				startpt[i] = 10.0*erand48(randBuffer)-5.0;
			}

			jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax, &funevals);
            #if DEBUG
                printf("\n\n\nHOOKE %d USED %d ITERATIONS, AND RETURNED\n", trial, jj);
                for (i = 0; i < nvars; i++)
                    printf("x[%3d] = %15.7le \n", i, endpt[i]);
            #endif

            //function evaluation
			fx = f(endpt, nvars, &funevals);
            #if DEBUG
                printf("f(x) = %15.7le\n", fx);
            #endif

			if (fx < best_fx) {
				best_trial = trial;
				best_jj = jj;
				best_fx = fx;
				for (i = 0; i < nvars; i++)
					best_pt[i] = endpt[i];
			}
		
		}

		//Initialize the input and output
		__attribute__((aligned(64))) double local_data[2] =  {funevals,best_fx};
		__attribute__((aligned(64))) double global_data[2] = {0,0};

		//create the operation
		MPI_Op MPI_MY_MIN;
		MPI_Op_create((MPI_User_function *)cust_min, 0, &MPI_MY_MIN);

		//Perform the reduction
		MPI_Allreduce(&local_data, &global_data, 2, MPI_DOUBLE, MPI_MY_MIN, MPI_COMM_WORLD);
			
		MPI_Op_free(&MPI_MY_MIN); //free the operator
		t1 = MPI_Wtime();
	
		if ( best_fx == global_data[1] ) {

			//Printing output to file
			FILE* output;
			char* stringOutput;
			if ( 0 > asprintf(&stringOutput,"MPI.%d.txt", size )) perror("String formatting failed"), MPI_Abort(MPI_COMM_WORLD,1);
			if ( (output = fopen(stringOutput , "w")) == NULL ) perror("Error at Accessing output file "), MPI_Abort(MPI_COMM_WORLD,1);

			//printing the results
			fprintf(output,"FINAL RESULTS: %d:%d\n",nvars,ntrials);
			fprintf(output,"Elapsed time = %.3lf s\n", t1-t0); //effective time 
			fprintf(output,"Total number of trials = %d\n", ntrials); //number of trials --just ntrials 
			fprintf(output,"Total number of function evaluations = %ld\n", (unsigned long)global_data[0]); //total function evaluations
			fprintf(output,"Total number of processes used = %d\n", size);
			fprintf(output,"Best result at trial %d of rank %d used %d iterations, and returned\n", best_trial, rank, best_jj);
			for (i = 0; i < nvars; i++) 
				fprintf(output,"x[%3d] = %15.7le \n", i, best_pt[i]);
			fprintf(output,"f(x) = %15.7le\n", best_fx);
		
			//free(output); 		//freeing the pointers here may lead to aborting writing to the file.
			//free(stringOutput); 	// I dont think that we encounter any memory leak. The processes end some steps further
		}
	MPI_Finalize();
	return 0;
}

