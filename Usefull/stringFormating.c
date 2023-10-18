#define __STDC_WANT_LIB_EXT2__ 1  //Define you want TR 24731-2:2010 extensions -- for gcc 

#include <stdio.h>
#include <stdlib.h>

int main(int argc , char ** argv){
    char* string; 
    int i=7;
    if(0> asprintf(&string,"Inserting a number into a string : %d\n", i)) perror("String formatting failed");
    printf("%s",string);

    free(string);
    return 0;
}