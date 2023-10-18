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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>


#define MAXVARS		(250)	/* max # of variables	     */
#define RHO_BEGIN	(0.9)	/* stepsize geometric shrink */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */

/* global variables */
unsigned long funevals = 0; //function evaluations (how many times function f was called)


/* Rosenbrock classic parabolic valley ("banana") function */
double f(double *x, int n)
{
    double fv;
    int i;

    funevals++;
    fv = 0.0;
    for (i=0; i<n-1; i++)   /* rosenbrock */
        fv = fv + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);

    return fv;
}

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double delta[MAXVARS], double point[MAXVARS], double prevbest, int nvars)
{
	double z[MAXVARS];
	double minf, ftmp;
	int i;
	minf = prevbest;				//initialize min f value as previous (=current) best value = fbefore (newf is being calculated in this function)
	for (i = 0; i < nvars; i++)
		z[i] = point[i];			//initialize z coordinates as newx's coordinates
	for (i = 0; i < nvars; i++) {
		z[i] = point[i] + delta[i]; //move point z to (newx + delta) point (delta was declared as startpoint[i] + rho (=step) or just rho)
		ftmp = f(z, nvars);			//calculate rosenbrock value on z
		if (ftmp < minf)			//if new value is lesser than minf = prevbest = fbefore, change minf to ftmp and procceed to next coordinate
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];		//change direction 
			z[i] = point[i] + delta[i];		//move point z
			ftmp = f(z, nvars);				//calculate rosenbrock value 
			if (ftmp < minf)
				minf = ftmp;				//if new value is lesser than minf = prevbest = fbefore, change minf to ftmp and procceed to next coordinate
			else
				z[i] = point[i];			//else return z to where it was (newx)
		}
	}
	for (i = 0; i < nvars; i++)				//in the end, no matter what has happened, newx is becoming the same as z. These changes leave the functions
		point[i] = z[i];

	return (minf);
}


int hooke(int nvars, double startpt[MAXVARS], double endpt[MAXVARS], double rho, double epsilon, int itermax)
{
	double delta[MAXVARS];
	double newf, fbefore, steplength, tmp;
	double xbefore[MAXVARS], newx[MAXVARS];
	int i, j, keep;
	int iters, iadj;

	for (i = 0; i < nvars; i++) {
		newx[i] = xbefore[i] = startpt[i];  //initialize newx and xbefore as the current starting point
		delta[i] = fabs(startpt[i] * rho);  //initialize delta coordinates with the absolute value 
		if (delta[i] == 0.0)				//of starting point's coordinates * step or just step (=rho)
			delta[i] = rho;  				
	}
	iadj = 0;
	steplength = rho;
	iters = 0;
	fbefore = f(newx, nvars);	//calculate rosenbrock function for newx (=current starting point)
	newf = fbefore;   			//store the value of rosenbrock function on newx (=current starting point), in variable newf
	while ((iters < itermax) && (steplength > epsilon)) {
		iters++;
		iadj++;
#if DEBUG
		printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);
		for (j = 0; j < nvars; j++)
			printf("   x[%2d] = %.4le\n", j, xbefore[j]);
#endif
		/* find best new point, one coord at a time */
		for (i = 0; i < nvars; i++) {
			newx[i] = xbefore[i];						   //initialize next point of search as the previous one (xbefore changes later in each iteration)
		}
		newf = best_nearby(delta, newx, fbefore, nvars);   //newx is also changed in there?
		/* if we made some improvements, pursue that direction */
		keep = 1;											//flag for no roundoff errors?
		while ((newf < fbefore) && (keep == 1)) {
			iadj = 0;										//no idea what it does
			for (i = 0; i < nvars; i++) {					//for each coordinate:
				/* firstly, arrange the sign of delta[] */
				if (newx[i] <= xbefore[i])
					delta[i] = 0.0 - fabs(delta[i]);		//if the  i-th coordinate of the new point is smaller than the pevious point, 
				else										//change direction on that coordinate
					delta[i] = fabs(delta[i]);				//else continus on the same direction on that coordinate
				/* now, move further in this direction */
				tmp = xbefore[i];							//store the i-th coordinate of xbefore (previous point of search) in temp
				xbefore[i] = newx[i];						//new point of search becomes the previous 
				newx[i] = newx[i] + newx[i] - tmp;			//i-th coordinate of newx is 2*(i-th coordinate of newx (= xbefore?) ) - (i-th coordinate of xbefore)
			}												                             
			fbefore = newf;									// newf is now fbefore
			newf = best_nearby(delta, newx, fbefore, nvars); 
			/* if the further (optimistic) move was bad.... */
			if (newf >= fbefore)
				break;															//break the inner while loop

			/* make sure that the differences between the new */
			/* and the old points are due to actual */
			/* displacements; beware of roundoff errors that */
			/* might cause newf < fbefore */

			keep = 0;															//assume no improvements were made
			for (i = 0; i < nvars; i++) {										//for each coordinate:
				keep = 1;														//assume improvement was made 
				if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i])))		//if |newx[i] - xbefore[i]| > 0.5*|delta[i]| , break the inner while loop (this expression probably means ???)
					break;														//break th inner while loop --> keep is still 1 - stop follwing the delta direction and go yo line 261
				else
					keep = 0;													//roundoff error occured --> don't pursue the same direction (in line 229)
			}
		}
		if ((steplength >= epsilon) && (newf >= fbefore)) {						//newf and fbefore from lines 242-243
			steplength = steplength * rho;										//What stands is rho = 0.9 < 1 (line 142), so steplength becomes smaller with each iteration 
			for (i = 0; i < nvars; i++) {
				delta[i] *= rho;												//change delta coordinates
			}
		}
	}
	for (i = 0; i < nvars; i++)			//initialize endpt point
		endpt[i] = xbefore[i];			//from line 239. These changes leave the function

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
	double startpt[MAXVARS], endpt[MAXVARS];
	int itermax = IMAX;
	double rho = RHO_BEGIN;
	double epsilon = EPSMIN;
	int nvars;
	int trial, ntrials;
	double fx;
	int i, jj;
	double t0, t1;

	double best_fx = 1e10;  //initialize absolute minimum
	double best_pt[MAXVARS];	//point on which the absolute minimum is found
	int best_trial = -1;	    //trial on which the absolute minimum is found
	int best_jj = -1;           //iterations on this trial

	for (i = 0; i < MAXVARS; i++) best_pt[i] = 0.0;	//initialize best point coordinates at [0,0,0 ... ,0]

	ntrials = 64*1024;	/* number of trials */ //= number of starting points (must be 64K = 2^16)
	nvars = 32;		/* number of variables (problem dimension) */
	srand48(1); //The initializer function srand48() sets the high-order 32 bits of Xi to the 32 bits contained in its argument. 
				//The low-order 16 bits of Xi are set to the arbitrary value 330E16 

	if (argc > 1) //set the ntrials as the first execution argument
    	{
       		ntrials = atoi(argv[1]);
    	}

	//write important results to file
	FILE* output;
	if ( (output = fopen( "sequential.txt", "w")) == NULL ) perror("Error at Accessing output file "), exit(1);

	t0 = get_wtime();
	for (trial = 0; trial < ntrials; trial++) {		 //amount of times a starting point will be created and the procedure will be followed4
		/* starting guess for rosenbrock test function, search space in [-5, 5) */
		for (i = 0; i < nvars; i++) {
			startpt[i] = 10.0*drand48()-5.0;		//starting point of search is selected (its coordinates are initialized) randomly 
		}

		//find ending point and store interations in jj
		jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);
#if DEBUG
		printf("\n\n\nHOOKE %d USED %d ITERATIONS, AND RETURNED\n", trial, jj);
		for (i = 0; i < nvars; i++)
			printf("x[%3d] = %15.7le \n", i, endpt[i]);
#endif

		fx = f(endpt, nvars); //calculate rosenbrock function's value on ending point
#if DEBUG
		printf("f(x) = %15.7le\n", fx);
#endif
		if (fx < best_fx) {
			best_trial = trial;  //initialize trial that gave the absolute minimum
			best_jj = jj;		 //initialize number of iterations on that trial
			best_fx = fx;		 //initialize absolute minimum value
			for (i = 0; i < nvars; i++)
				best_pt[i] = endpt[i];  //initialize point on which the absolute minimum was found
		}
	}
	
	t1 = get_wtime();
		
	fprintf(output,"FINAL RESULTS: %d:%d\n",nvars,ntrials);
	fprintf(output,"Elapsed time = %.3lf s\n", t1-t0); //effective time 
	fprintf(output,"Total number of trials = %d\n", ntrials); //number of trials --just ntrials 
	fprintf(output,"Total number of function evaluations = %ld\n", funevals); //total function evaluations
	fprintf(output,"Best result at trial %d used %d iterations, and returned\n", best_trial, best_jj);
	for (i = 0; i < nvars; i++) {
		fprintf(output,"x[%3d] = %15.7le \n", i, best_pt[i]);
	}
	fprintf(output,"f(x) = %15.7le\n", best_fx);

	//free(output); //leads errors
    return 0;
}
