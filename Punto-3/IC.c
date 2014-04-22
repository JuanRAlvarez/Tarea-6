#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 4.296E-6 //units of (km*2 * kpc) msun*-1 s*-2

FLOAT *get_memory(int n_points);

int main(int argc, char **argv){
    
    
    //Obtengo los parámetros del mundo real. 1,2,3 posiciones x,y,z. 4,5,6 velocidades iniciales 7,8,9 semilla, masa punto central, radio círculo, 10 es el número de partículas
    
    int i;
    
    FLOAT x0 = atof(argv[1]);
    FLOAT y0 = atof(argv[2]);
    FLOAT z0 = atof(argv[3]);
    
    FLOAT v_0x = atof(argv[4]);
    FLOAT v_0y = atof(argv[5]);
    FLOAT v_0z = atof(argv[6]);
    
    int seed = atoi(argv[7]);
    FLOAT M = atof(argv[8]);
    FLOAT R = atof(argv[9]);
    int N = atoi(argv[10]);
    
    //
    
    srand48(seed);
    
    FLOAT x;
    FLOAT y;
    FLOAT z;
    FLOAT r;
    FLOAT vel;
    
    FLOAT v_x;
    FLOAT v_y;
    FLOAT v_z;
    FLOAT mass;
    
    int ID = -1;
    
    FILE *in;
    
    FLOAT erre;
    FLOAT theta;
    
    in = fopen("IC.dat","w");
    
    fprintf(in, "%i %f %f %f %f %f %f %f \n", ID, x0, y0, z0, v_0x,v_0y,v_0z, M);
    
    for (i=0; i<N; i++) {
        
        ID = i+1;
        erre = drand48() * R;
        theta = 2*PI*drand48();
        
        
        //x = pow(erre,0.5)*cos(theta) + x0;
        //y = pow(erre,0.5)*sin(theta) + y0;
        x = cos(theta) + x0;
        y = sin(theta) + y0;
        z = 0 + z0;
        
        r = sqrt(pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2));
        
        if (r!=0) {
            
            
            vel = sqrt(G_GRAV * M / r);
            
            v_x = v_0x - vel*(y-y0)/r;
            v_y = v_0y + vel*(x-x0)/r;
            v_z = v_0z;
        }
        else{
            v_x = v_0x - vel;
            v_y = v_0y;
            v_z= v_0z;
        }
        
        
        
        mass = 1.0;
        
         fprintf(in, "%i %f %f %f %f %f %f %f \n", ID, x, y, z, v_x,v_y,v_z, mass);
        
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