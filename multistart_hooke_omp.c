/* Nonlinear Optimization using the algorithm of Hooke and Jeeves  */
/*	12 February 1994	author: Mark G. Johnson 	   */


/* Find a point X where the nonlinear function f(X) has a local    */
/* minimum.  X is an n-vector and f(X) is a scalar.  In mathe-	   */
/* matical notation  f: R^n -> R^1.  The objective function f()    */
/* is not required to be continuous.  Nor does f() need to be	   */
/* differentiable.  The program does not use or require 	   */
/* derivatives of f().						   */

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
/*  FITNESS FOR ANY PARTICULAR PURPOSE. 			   */
/*								   */

#define __STDC_WANT_LIB_EXT2__ 1  //Define you want TR 24731-2:2010 extensions -- for gcc 

#include <stdio.h> // I/O operations
#include <stdlib.h> //memory allocation and randomizers
#include <math.h> //needs -lm at the end
#include <time.h> // Time struct for open mp
#include <omp.h> //include openmp
#include <sys/time.h> // Time struct
#include <string.h> //memset

unsigned long funevals = 0; // number of function evaluations (how many times function f was called)
							// Total iterations count = total calls of rosenbrocks parabolic valley
							// Left global for the debugging purposes

//set openmp environmental variables -- IN the terminal
//export OMP_DYNAMIC=FALSE -- or omp_set_dynamic(0) = false and omp_set_dynamic(1) = true
//export OMP_PROC_BIND=TRUE

#define MAXVARS		(32)	/* max # of variables	     */  //The maximum number of variables 
#define RHO_BEGIN	(0.9)	/* stepsize geometric shrink */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */ // number of trials 

/* Rosenbrock classic parabolic valley ("banana") function */
double f(double *x, int n, unsigned long *l_funevals) //endpt , and nvars 
{
    double fv = 0;
    int i;

    (*l_funevals)++; 	// this not the  augmentation of the global counter of function evaluations but it is localized
						// for each thread
    // n is nvars
    for (i=0; i<n-1; i++)   /* rosenbrock */ //just 31 iterations
        fv = fv + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);

    return fv;
}

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double delta[MAXVARS], double point[MAXVARS], double prevbest, int nvars, unsigned long* l_funevals)
{
	double z[MAXVARS];
	double minf, ftmp;
	int i;
	minf = prevbest;
	for (i = 0; i < nvars; i++) 
		z[i] = point[i];
	for (i = 0; i < nvars; i++) {
		z[i] = point[i] + delta[i];
		ftmp = f(z, nvars,l_funevals);
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i] = point[i] + delta[i];
			ftmp = f(z, nvars,l_funevals);
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

//nvars is the number of variables == problem size

int hooke(int nvars, double startpt[MAXVARS], double endpt[MAXVARS], double rho, double epsilon, int itermax,unsigned long *l_funevals)
{
	double delta[MAXVARS];
	double newf, fbefore, steplength, tmp;
	double xbefore[MAXVARS], newx[MAXVARS];
	int i, j, keep; //j is used for debug mode
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
	fbefore = f(newx, nvars, l_funevals);
	newf = fbefore;
	while ((iters < itermax) && (steplength > epsilon)) {
		iters++;
		iadj++;
#if DEBUG
		printf("\nAfter %5lu funevals, f(x) =  %.4le at\n", funevals, fbefore);
		for (j = 0; j < nvars; j++)
			printf("   x[%2d] = %.4le\n", j, xbefore[j]);
#endif
		/* find best new point, one coord at a time */
		for (i = 0; i < nvars; i++) {
			newx[i] = xbefore[i];
		}
		newf = best_nearby(delta, newx, fbefore, nvars, l_funevals);
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
			newf = best_nearby(delta, newx, fbefore, nvars, l_funevals);
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



double get_wtime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}

int main(int argc, char *argv[])
{
	int itermax = IMAX; 		//read only initialized variables
	double rho = RHO_BEGIN; 	//
	double epsilon = EPSMIN; 	//
	int nvars = 32; 			//number of variables 
	int ntrials = 64*1024;		//number of trials (problem dimentions)
	
	double t0, t1;				//timers
	int k, threads=1; 			//general counter for main, and default threads used

	double best_fx = 1e10;		//reduction		//initialize absolute minimum
	double best_pt[MAXVARS];	//reducton		//point on which the absolute minimum is found
	//Initializing the Array is unecessary and justs adds to our time complexitiy

	int best_trial = -1;		//initialization of best τrial where min is found
	int best_jj = -1;			//initializaton of the counter of iterations for the best trial
	long l_timeSeed = time(0);

	if (argc > 2) //set ntrials and thread number as the first and second execution argument respectively
    {
		ntrials = atoi(argv[1]);
	   	threads = atoi(argv[2]);
    }

	//OMP ENVIRONMENTAL VARIABLES
	omp_set_num_threads(threads); 	//initialize thread number as 1 if no argument has been given
	omp_set_dynamic(0); 			//We want the explicit amount of threads we are requesting

	//writing to file instead of console 
	FILE* output;
	char* stringOutput;
	if ( 0> asprintf(&stringOutput,"Omp.%d.txt", threads)) 	perror("String formatting failed"), exit(1);
	if ( (output = fopen(stringOutput , "w")) == NULL ) 	perror("Error at Accessing output file "), exit(1);

	t0 = get_wtime(); //starting the clock

 	/*------	Here we start the main paralization 	------*/
	#pragma omp parallel	//argv1 num threads 
	{
		//Declared variables in this scope are private for each thread

		__attribute__((aligned(64))) double startpt[MAXVARS];  
		__attribute__((aligned(64))) double endpt[MAXVARS];

		int i, jj, trial;				//i is a counter , jj the result of hooke and trial signifies the current trial of each thread
		double fx;						//temporary rosenbrock fuction evaluation for each loop
		unsigned long l_funevals = 0; 	//local funection evaluations of each thread

		//Initializing local best iterations
		double l_best_fx = 1e10;
		double l_best_pt[MAXVARS];
		double l_best_jj = -1;
		double l_best_trial= -1;
		//End of moved declarations 

		unsigned short randBuffer[3] = {0,0,l_timeSeed+omp_get_thread_num()}; //just 6 bytes no need for alignement

		//We opt for static scheduling for the better memory managerment. We do not deem that the workload of each loop changes as the number of iterations is augmentated.
		//The difference in the workload of each loop is actually determined by the pseudo-random initialization of the startpt array.
		//And as it is pseudo-random we are guessing that changing the scheduling does not add any concrete advantage.

		#pragma omp for schedule(static) nowait		//nowait is used to skip the implicit barrier at the end of <omp for> 
		for (trial = 0; trial < ntrials; trial++)   //amount of times a starting point will be created and the procedure will be followed
		{
			//first pseudo-random guess for rosenbrock test function, search space in [-5, 5) 
			for (i = 0; i < nvars; i++) { 
				//erand48() gives values of (0,1)
				startpt[i] = 10.0*erand48(randBuffer)-5.0;
			}
						
			//find ending point and store it in jj
			jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax, &l_funevals);
			#if DEBUG
				printf("\n\n\nHOOKE %d USED %d ITERATIONS, AND RETURNED\n", trial, jj);
				for (i = 0; i < nvars; i++)
					printf("x[%3d] = %15.7le \n", i, endpt[i]);
			#endif

			//calculate rosenbrock function's value on ending point
			fx = f(endpt, nvars, &l_funevals); 

			#if DEBUG
				printf("f(x) = %15.7le\n", fx);
			#endif

			//best trial of each thread
			if (fx < l_best_fx) {
				//store the data of best iteration 
				l_best_trial = trial; 
				l_best_jj = jj;	
				l_best_fx = fx;	
				for (i = 0; i < nvars; i++) 
					l_best_pt[i] = endpt[i];
			}
		}

		// Reduction of best trial between threads
		//creates a critical section to avoid race conditions between the threads
		#pragma omp critical 
		{
			if (l_best_fx < best_fx) {
				//store the absolute best data of each thread 
				best_trial = l_best_trial;
				best_jj = l_best_jj;	
				best_fx = l_best_fx;	
				for (i = 0; i < nvars; i++) 
					best_pt[i] = l_best_pt[i]; 
			}
			funevals += l_funevals; 
		}

	} //end of omp parallel 
 	/*------	Here we end the main paralization ------*/

	t1 = get_wtime(); 
	
	//writing the results to output file
	fprintf(output,"FINAL RESULTS: %d:%d\n",nvars,ntrials);
	fprintf(output,"Elapsed time = %.3lf s\n", t1-t0); 							//effective time 
	fprintf(output,"Total number of trials = %d\n", ntrials); 					//ntrials 
	fprintf(output,"Total number of function evaluations = %ld\n", funevals); 	//total function evaluations
	fprintf(output,"Total number of threads used = %d\n", threads); 			//threads used
	fprintf(output,"Best result at trial %d used %d iterations, and returned\n", best_trial, best_jj);
	for (k = 0; k < nvars; k++) {
		fprintf(output,"x[%3d] = %15.7le \n", k, best_pt[k]);
	}
	fprintf(output,"f(x) = %15.7le\n", best_fx);

	//free(output); //some times the results are not written before the freeing of the FILE*
	free(stringOutput);
	return 0;
}
