#!/bin/bash

# The script is used to run the makefile and the produced executables once

#Execute script using 3 arguments (unsigned integers ): trialsNum threadNum OptimizationLevel
#trialsNum  -> ntrials (default 64k)
#ThreadNum  -> Threads or Processes that will be used
#OptimizationLevel  -> Compilation options :
#   Optimization 1 -> -O3 
#   Optimization 2 -> -O3 -mfma -msse2
#   Optimization 3 -> -O3 -march=native -flto -ffast-math -fprofile-generate

#make clean -> remove executables
#make cleantxt -> remove results

function runExecutables {
    export OMP_PROC_BIND=$bool && ./mult_hooke_omp $trials $threads
    export OMP_PROC_BIND=$bool && ./mult_hooke_tasks $trials $threads
    mpiexec -n $threads --oversubscribe mult_hooke_mpi $trials 

    if [ $threads -eq 1 ]; then
        ./mult_hooke_seq $trials
        export OMP_PROC_BIND=$bool && mpiexec -n $threads --oversubscribe mult_hooke_mpi_omp $trials 1
    else
        hthreads=$(($threads/2))
        export OMP_PROC_BIND=$bool && mpiexec -n $hthreads --oversubscribe mult_hooke_mpi_omp $trials $hybridThreads
    fi
}
#Change the number of the Threads used by each process in the hybrid ( MPI + OPENMP ) parallelization 
hybridThreads=2
bool=false
make cleanall

if [ $# -ne 3 ]; then
    echo "Execute script using 3 arguments (unsigned integers ): trialsNum MaxThreadNum OptimizationLevel -> ( 1 or 2 or 3)"
    echo "Execution follows default values: 64k trials, 14 MaxThreadNum, 0 OptimizationLevel"

    #Initialize execution arguments
    trials=64000
    threads=6
    make
    runExecutables
    echo "Scipt Execution time = $SECONDS s"
 
else
    #Initialize execution arguments
    trials=$1
    threads=$2
    make aggresive_optimization=$3
    runExecutables
    echo "Scipt Execution time = $SECONDS s"

fi

#Comments about Bash scripting :
# Use $# for the number of execution arguments
# Access to the execution variables is through $1 , $2 etc...
