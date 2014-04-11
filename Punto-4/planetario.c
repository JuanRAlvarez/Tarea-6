#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define G 39.486 //(au)^3 /(m(+) y)
#define PI 3.14159265

/* ESTE CÓDIGO ESTÁ BASADO EN EL CÓDIGO DEL PROFESOR JAIME FORERO (github.com/forero), PERO TIENE ALGUNAS MODIFICACIONES, A SABER:
 se usan sólo dos dimensiones
 se usa una sola función para inicializarlo todo
 se implementa el método de runge kutta de cuarto orden
 */


FLOAT *get_memory(int n_points);
void initialize_everything(FLOAT *x, FLOAT *y, FLOAT *vx, FLOAT *vy, FLOAT *ax, FLOAT *ay, int n_points, FLOAT vel, FLOAT radius,FLOAT * mass, FLOAT unit_mass);
void print_velocities(FLOAT *x, FLOAT *y, FLOAT *vx, FLOAT *vy, FLOAT *ax, FLOAT *ay, int n_points);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *vx, FLOAT *vy, FLOAT *x,FLOAT *y, FLOAT *mass, int n_points);
void escribir(float *x, float *y, int pasos, char *filename);

/*
 Unidades:
 Masa: masas solares
 Distancia: Unidades atronómicas
 Tiempo: años
 */

int main(int argc,char **argv){

    int i,j;
    
    /*positions of all particles*/
    FLOAT *x;
    FLOAT *y;
    
    /*velocities of all particles*/
    FLOAT *vx;
    FLOAT *vy;
    
    /*accelerations of all particles*/
    FLOAT *ax;
    FLOAT *ay;
    
    /*masses*/
    FLOAT *mass;
    FLOAT unit_mass = 1.0;
    
    /*Timestep variables*/
    FLOAT delta_t = 0.1;
    int n_points = 3;
    int n_steps = (int)(5000.0/delta_t);
    FLOAT radius = 100.0;
    
    /*memory allocation*/
    x = get_memory(n_points);
    y = get_memory(n_points);
    vx = get_memory(n_points);
    vy = get_memory(n_points);
    ax = get_memory(n_points);
    ay = get_memory(n_points);
    mass = get_memory(n_points);
    
    FLOAT vel = sqrt((11.0/3.0) * G * unit_mass / (sqrt(3.0)*radius));
    
    initialize_everything(x, y, vx, vy, ax, ay, n_points, vel, radius, mass, unit_mass);
    
    /*implementation of a simple Euler integration*/
    for(i=0;i<n_steps;i++){
        get_acceleration(ax, ay,vx,vy, x, y, mass, n_points);
        for(j=0;j<n_points;j++){
            x[j] = x[j] + delta_t * vx[j];
            y[j] = y[j] + delta_t * vy[j];
            
            vx[j] = vx[j] + delta_t * ax[j];
            vy[j] = vy[j] + delta_t * ay[j];
        }
        print_velocities(x, y, vx, vy, ax, ay, n_points);
    }
    
    return 0;
}


void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *vx, FLOAT *vy, FLOAT *x,FLOAT *y, FLOAT *mass, int n_points){
    
    int i,j;
    FLOAT r_ij;
    
    for(i=0;i<n_points;i++){
        ax[i]=0.0;
        ay[i]=0.0;
        
        for(j=0;j<n_points;j++){
            if(j!=i){
                r_ij = (pow((x[i] - x[j]),2.0) +
                        pow((y[i] - y[j]),2.0) );
                r_ij = sqrt(r_ij);
                ax[i] += -G *mass[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
                ay[i] += -G *mass[j]/ pow(r_ij,1.5) * (y[i] - y[j]);
            }
        }
    }
}


FLOAT * get_memory(int n_points){
    FLOAT * x;
    if(!(x = malloc(sizeof(FLOAT) * n_points))){
        printf("problem with memory allocation");
        exit(1);
    }
    return x;
}

void initialize_everything(FLOAT *x, FLOAT *y, FLOAT *vx, FLOAT *vy, FLOAT *ax, FLOAT *ay, int n_points, FLOAT vel, FLOAT radius,FLOAT * mass, FLOAT unit_mass){

    int i;

    FLOAT angle;
    angle = 2.0*PI/n_points;
    for (i=0; i<n_points; i++) {
    
    x[i]= cos(angle * i) * radius;
    y[i]= sin(angle * i )* radius;
    
    vx[i] = -sin(angle * i) * vel;
    vy[i] = cos(angle * i) * vel;
        
    ax[i] = 0.0;
    ay[i] = 0.0;
        
    }
    
}


void print_velocities(FLOAT *x, FLOAT *y, FLOAT *vx, FLOAT *vy, FLOAT *ax, FLOAT *ay, int n_points){

    int i;
    for(i=0;i<n_points;i++){
        printf("%f %f \n", vx[i], vy[i]);
    }

}

void escribir(float *x, float *y, int pasos, char *filename){
    
    
    FILE *in;
    
    in = fopen(filename,"w");
    int i;
    
    
    fprintf(in,"%f %f %f %f %f %f \n", x[1], y[1], x[2], y[2], x[3], y[3]);
    
    fclose(in);
    
}
