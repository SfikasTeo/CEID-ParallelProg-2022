#define _GNU_SOURCE  //The macro is necessary

#include <omp.h>
#include <stdio.h>
#include <sched.h>


void test_proc_bind()
{
#pragma omp parallel
    {
    int tid = omp_get_thread_num();
    int core =  sched_getcpu(); /* linux specific */
    #pragma omp critical
        printf("Thread %d running on core %d\n", tid, core);
    }
}


int main()
{
    omp_set_num_threads(4);
    test_proc_bind();
    return 0;
}
