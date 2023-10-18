#include <sys/time.h>
#include <stdio.h>
//#include <windows.h> // for core windows utilities like Sleep
#include <unistd.h> // for linux

// timer
double get_wtime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1.0e-6;
}

int main(int argc , char ** argv){
    
    const double t0 = get_wtime();
    //Sleep(3000); //3000 milliseconds for windows
    sleep(3); // 3 seconds for linux
    const double t1 = get_wtime();

    const double teval = t1 - t0;
    printf("Execution time is : %f", teval);
}