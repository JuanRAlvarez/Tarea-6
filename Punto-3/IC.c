#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1

FLOAT *get_memory(int n_points);

int main(int argc, char **argv){
    
    srand48(100);
    
    int i;
    
    FLOAT x0 = atof(argv[1]);
    FLOAT y0 = atof(argv[2]);
    FLOAT z0 = atof(argv[3]);
    
    FLOAT v_0x = atof(argv[4]);
    FLOAT v_0y = atof(argv[5]);
    FLOAT v_0z = atof(argv[6]);
    
    FLOAT M = atof(argv[7]);
    FLOAT R = atof(argv[8]);
    int N = atoi(argv[9]);
    
    FLOAT x;
    FLOAT y;
    FLOAT z;
    
    FLOAT v_x;
    FLOAT v_y;
    FLOAT v_z;
    
    int ID = -1;
    
    FILE *in;
    
    in = fopen("IC.dat","w");
    
    fprintf(in, "%i %f %f %f %f %f %f \n", ID, x0, y0, z0, v_0x,v_0y,v_0z);
    
    ID = 1;
    for (i=0; i<N; i++) {
        x = R*2*(drand48()-0.5) + x0;
        y = R*2*(drand48()-0.5) + y0;
        z = 0 + z0;
        
        v_x = -y;
        v_y = x;
        v_z = 0;
        
         fprintf(in, "%i %f %f %f %f %f %f \n", ID, x, y, z, v_x,v_y,v_z);
        
    }
    
    
    fclose(in);
    
    
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