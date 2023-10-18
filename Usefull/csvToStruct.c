#include <stdio.h> // default I/O operations 
#include <stdlib.h> // memory management
#include <string.h>

#define MAX_STR_LEN 256 // incomplete and dangerous memory management
#define MAX_RECORD_NUM 5 //statical could be done dynamically -> see below 

typedef struct { // records from csv
    char* username;
    char* identifier;
    char* firstName;
    char* lastName; 
}record;

/*The Dataframe is a double pointer to record for simplicity reasons -- void ** should be used for generalized applications*/
record** readCSVtoDataframe(char* csv){ // on call may argv[1] be used
    // -- Memory Allocation -- 
    record ** dataframe; // Array of record pointers
    if ( (dataframe = (record**)malloc(sizeof(record*) * MAX_RECORD_NUM)) == NULL ) perror("Insufficient Memory"), exit(EXIT_FAILURE);// memory allocation
    char* tmp;
    if ( (tmp = malloc(MAX_STR_LEN)) == NULL) perror("Insufficient Memory"), exit(EXIT_FAILURE); 

    // -- Open the I/O stream
    FILE* buffer;
    if ( (buffer = fopen(csv,"rb")) == NULL ) perror("Error at opening CSV"), exit(EXIT_FAILURE);

    /* Initialization of the dataframe from the csv */
    
    char* nextcolumn;
    for (int i=0; i<=MAX_RECORD_NUM; i++){
        fgets(tmp,255,buffer);
        if((strlen(tmp)>0) && (tmp[strlen(tmp)-1] == '\n')) tmp[strlen(tmp)-1] = '\0'; // if the length of the line is greater than one and ends with newline character -> make last character EOF of the tmp string
        if((*(dataframe+i) = (record*)malloc(sizeof(record))) == NULL) perror("Insufficient Memory"), exit(EXIT_FAILURE);  // manage memory allocation for each structure
    
        nextcolumn = strtok(tmp, ";"); //delimiter of the CSV file is ";"
        (*(*(dataframe+i))).username = strdup(nextcolumn); //We use strdup to create the needded allocated space for the dataframe member variable (string)
        /* dataframe+i -> memory address of a record || *(dataframe+i) -> pointer of a record ||  (*(*(dataframe+i))).username -> (the record i).username */

        nextcolumn = strtok(NULL, ";");
        (*(*(dataframe+i))).identifier = strdup(nextcolumn);
        nextcolumn = strtok(NULL, ";");
        (*(*(dataframe+i))).firstName = strdup(nextcolumn);
        nextcolumn = strtok(NULL, ";");
        (*(*(dataframe+i))).lastName = strdup(nextcolumn);
        
    }
    fclose(buffer);
return dataframe;
}

int main(int argc , char** argv){
    if (argc > 1) { 
        record ** dataframe = readCSVtoDataframe(argv[1]);
        for (int i = 0; i<= MAX_RECORD_NUM; i++ )
        printf("index i= %i  ID: %s, %s, %s, %s \n",i, dataframe[i]->username , dataframe[i]->identifier, dataframe[i]->firstName , dataframe[i]->lastName);
    }
    
return 0;
}

// Warnings -- memory leak may exist but i dont think so. 
// 1) making a array of pointers to structs may be helpfull but not always -- fun to make but we ought to be carefull of cache misses.
// For most of the data sizes that fit easily into to cashe ( and especially for small structures ) using struct[] -> struct array will be easier and more efficient.
// Wrong but a helpfull implementation follows :  https://stackoverflow.com/questions/20212714/reading-a-csv-file-into-struct-array
// Better but more advanced solution : https://ideone.com/mSCgPM (no struct)
// 2) Using fread as follows : 
    /* 
        fseek(file , 0L , SEEK_END); //moving file pointer to the end
        lsize = ftell(file); //returns current file position
        rewind(file); // resets file pointer

        // allocate heap memory for the input
        buffer = calloc(1,lsize+1);
        if( !buffer ) fclose(file),fputs("memory alloc fails",stderr),exit(1);

        // copy the file into the buffer 
        if( 1!=fread( buffer , lsize, 1 , file) )
         fclose(file),free(buffer),fputs("entire read fails",stderr),exit(1);
   
        //close the FILE 
        fclose(file);  
        //PROCESSING OF DATA
        free(buffer);
    */ 
// -> and processing the data through the buffered state may increase exponentially the speedup in constrast with the stream. Depends on the data size. 
// This also may help with making the MAX_REC_NUM dynamic