#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1

FLOAT *get_memory(int n_points);

int main(int argc, char **argv){
    
    return 0;
    
}

FLOAT *get_memory(int n_points){
    FLOAT * x;
    if(!(x = malloc(sizeof(FLOAT) * n_points))){
        printf("problem with memory allocation");
        exit(1);
    }
    return x;
}