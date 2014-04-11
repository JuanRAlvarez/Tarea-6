#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1

FLOAT * get_memory(int n_points);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points, int *ID);
FLOAT get_kinetic(FLOAT *v_x, FLOAT *v_y, FLOAT *v_z, FLOAT *mass, int n_points);
FLOAT get_potential(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);

int main(int argc, char **argv){
    
    //Obtengo los parámetros del mundo exterior
    FLOAT T = atof(argv[1]);
    int M = atoi(argv[2]);
    int n_points = atoi(argv[3]);
    
    
    /*accelerations of all particles*/
    FLOAT *a_x;
    FLOAT *a_y;
    FLOAT *a_z;
    
    /*terms for runge-kutta*/
    
    FLOAT *a_x_old;
    FLOAT *a_y_old;
    FLOAT *a_z_old;
    
    
    FLOAT  *k_1_x;
    FLOAT  *k_1_y;
    FLOAT  *k_1_z;
    
    FLOAT  *k_1_v_x;
    FLOAT  *k_1_v_y;
    FLOAT  *k_1_v_z;
    
    FLOAT  *k_2_x;
    FLOAT  *k_2_y;
    FLOAT  *k_2_z;
    
    FLOAT  *k_2_v_x;
    FLOAT  *k_2_v_y;
    FLOAT  *k_2_v_z;
    
    FLOAT  *k_3_x;
    FLOAT  *k_3_y;
    FLOAT  *k_3_z;
    
    FLOAT  *k_3_v_x;
    FLOAT  *k_3_v_y;
    FLOAT  *k_3_v_z;
    
    FLOAT  *k_4_x;
    FLOAT  *k_4_y;
    FLOAT  *k_4_z;
    
    FLOAT  *k_4_v_x;
    FLOAT  *k_4_v_y;
    FLOAT  *k_4_v_z;
    
    FLOAT *xtemp;
    FLOAT *ytemp;
    FLOAT *ztemp;
    FLOAT *v_xtemp;
    FLOAT *v_ytemp;
    FLOAT *v_ztemp;
    FLOAT *a_xtemp;
    FLOAT *a_ytemp;
    FLOAT *a_ztemp;
    
    /*NOTA! TODAVÍA NO SÉ CÓMO CALCULAR EL TAMAÑO DEL ARCHIVO QUE VIENE COMO PARÁMETRO DE ENTRADA. POR AHORA, ESTOY PENSANDO DEJARLO COMO UN PARÁMETRO DE ENTRADA DE ESTE ARCHIVO
     Y PROCEDER A INGRESAR ESTE VALOR CUANDO TENGA UN MAKEFILE*/
    
    int i;
    
    int *ID;
    FLOAT *x;
    FLOAT *y;
    FLOAT *z;
    FLOAT *v_x;
    FLOAT *v_y;
    FLOAT *v_z;
    
    /*masses*/
    FLOAT *mass;
    
    /*timestep variables*/
    FLOAT h= 0.001;
    int n_steps = (int)(100.0/h);
    int n_points = 3;
    FLOAT radius = 100.0;
    FLOAT unit_mass = 1.0;
    FLOAT vel_initial = sqrt((11.0/3.0) * G_GRAV * unit_mass / (sqrt(3.0)*radius));
    int i,j,k;
    
    /*memory allocation*/
    x = get_memory(n_points);
    y = get_memory(n_points);
    z = get_memory(n_points);
    v_x = get_memory(n_points);
    v_y = get_memory(n_points);
    v_z = get_memory(n_points);
    a_x = get_memory(n_points);
    a_y = get_memory(n_points);
    a_z = get_memory(n_points);
    mass = get_memory(n_points);
    a_x_old = get_memory(n_points);
    a_y_old = get_memory(n_points);
    a_z_old = get_memory(n_points);
    k_1_x = get_memory(n_points);
    k_1_y = get_memory(n_points);
    k_1_z = get_memory(n_points);
    k_1_v_x = get_memory(n_points);
    k_1_v_y = get_memory(n_points);
    k_1_v_z = get_memory(n_points);
    k_2_x = get_memory(n_points);
    k_2_y = get_memory(n_points);
    k_2_z = get_memory(n_points);
    k_2_v_x = get_memory(n_points);
    k_2_v_y = get_memory(n_points);
    k_2_v_z = get_memory(n_points);
    k_3_x = get_memory(n_points);
    k_3_y = get_memory(n_points);
    k_3_z = get_memory(n_points);
    k_3_v_x = get_memory(n_points);
    k_3_v_y = get_memory(n_points);
    k_3_v_z = get_memory(n_points);
    k_4_x = get_memory(n_points);
    k_4_y = get_memory(n_points);
    k_4_z = get_memory(n_points);
    k_4_v_x = get_memory(n_points);
    k_4_v_y = get_memory(n_points);
    k_4_v_z = get_memory(n_points);
    xtemp = get_memory(n_points);
    ytemp = get_memory(n_points);
    ztemp = get_memory(n_points);
    v_xtemp = get_memory(n_points);
    v_ytemp = get_memory(n_points);
    v_ztemp = get_memory(n_points);
    a_xtemp = get_memory(n_points);
    a_ytemp = get_memory(n_points);
    a_ztemp = get_memory(n_points);
    
    FILE *in1;
    in1 = fopen("IC.dat","r");
    
    //VOY A LEER EL ARCHIVO. POR AHORA ME SALE UN SEGMENTATION FAULT. :(
    
    /*for (i=0; i<n_points; i++) {
        fscanf(in1,"%i %f %f %f %f %f %f \n", &(ID[i]),&(x[i]),&(y[i]),&(z[i]),&(v_x[i]),&(v_y[i]),&(v_z[i]));
    }
    
    fclose(in1);
    
     
    
    for (i=0; i<10; i++) {
        printf("%i\n",ID[i]);
    }
    
    return 0;*/
    
    printf("%f %i %i \n",T,M, n_points);
    
    
    /*EL RUNGE-KUTTAZO ESTÁ INSPIRADO EN EL DEL MAIN DE 3 CUERPOS*
     
     /*implementation of a second order runge kutta integration*/
    
    FILE *in;
    in = fopen("3cuerpos.dat","w");
    
    FLOAT *x_old;
    FLOAT *y_old;
    
    for(i=0;i<n_steps;i++){
        
        a_x_old = a_x;
        a_y_old = a_y;
        a_z_old = a_z;
        
        get_acceleration(a_x, a_y, a_z, x, y, z, mass, n_points);
        for(j=0;j<n_points;j++){
            
            k_1_x[j] = v_x[j];
            k_1_y[j] = v_y[j];
            k_1_z[j] = v_z[j];
            
            k_1_v_x[j] = a_x[j];
            k_1_v_y[j] = a_y[j];
            k_1_v_z[j] = a_z[j];
            
            
            //FIRST STEP
            
            xtemp[j] = x[j] + (h/2.0)*k_1_x[j];
            ytemp[j] = y[j] + (h/2.0)*k_1_y[j];
            ztemp[j] = z[j] + (h/2.0)*k_1_z[j];
            
            v_xtemp[j] = v_x[j] + (h/2.0)*k_1_v_x[j];
            v_ytemp[j] = v_y[j] + (h/2.0)*k_1_v_y[j];
            v_ztemp[j] = v_z[j] + (h/2.0)*k_1_v_z[j];
            
            k_2_x[j] = v_xtemp[j];
            k_2_y[j] = v_ytemp[j];
            k_2_z[j] = v_ztemp[j];
            
        }
        
        get_acceleration(a_xtemp,a_ytemp,a_ztemp,xtemp,ytemp,ztemp, mass, n_points);
        
        for(j=0;j<n_points;j++){
            
            k_2_v_x[j] = a_xtemp[j];
            k_2_v_y[j] = a_ytemp[j];
            k_2_v_z[j] = a_ztemp[j];
            
            
            //SECOND STEP
            
            xtemp[j] = x[j] + (h/2.0)*k_2_x[j];
            ytemp[j] = y[j] + (h/2.0)*k_2_y[j];
            ztemp[j] = z[j] + (h/2.0)*k_2_z[j];
            
            v_xtemp[j] = v_x[j] + (h/2.0)*k_2_v_x[j];
            v_ytemp[j] = v_y[j] + (h/2.0)*k_2_v_y[j];
            v_ztemp[j] = v_z[j] + (h/2.0)*k_2_v_z[j];
            
            k_3_x[j] = v_xtemp[j];
            k_3_y[j] = v_ytemp[j];
            k_3_z[j] = v_ztemp[j];
            
        }
        
        get_acceleration(a_xtemp,a_ytemp,a_ztemp,xtemp,ytemp,ztemp, mass, n_points);
        
        for(j=0;j<n_points;j++){
            
            k_3_v_x[j] = a_xtemp[j];
            k_3_v_y[j] = a_ytemp[j];
            k_3_v_z[j] = a_ztemp[j];
            
            
            //THIRD STEP
            
            xtemp[j] = x[j] + h*k_3_x[j];
            ytemp[j] = y[j] + h*k_3_y[j];
            ztemp[j] = z[j] + h*k_3_z[j];
            
            v_xtemp[j] = v_x[j] + h*k_3_v_x[j];
            v_ytemp[j] = v_y[j] + h*k_3_v_y[j];
            v_ztemp[j] = v_z[j] + h*k_3_v_z[j];
            
            k_4_x[j] = v_xtemp[j];
            k_4_y[j] = v_ytemp[j];
            k_4_z[j] = v_ztemp[j];
            
        }
        
        get_acceleration(a_xtemp,a_ytemp,a_ztemp,xtemp,ytemp,ztemp, mass, n_points);
        
        for(j=0;j<n_points;j++){
            
            k_4_v_x[j] = a_xtemp[j];
            k_4_v_y[j] = a_ytemp[j];
            k_4_v_z[j] = a_ztemp[j];
            
            
            //FOURTH STEP
            
            x[j] = x[j] + h*(1.0/6.0)*(k_1_x[j]+2*k_2_x[j]+2*k_3_x[j]+k_4_x[j]);
            y[j] = y[j] + h*(1.0/6.0)*(k_1_y[j]+2*k_2_y[j]+2*k_3_y[j]+k_4_y[j]);
            z[j] = z[j] + h*(1.0/6.0)*(k_1_z[j]+2*k_2_z[j]+2*k_3_z[j]+k_4_z[j]);
            
            v_x[j] = v_x[j] + h*(1.0/6.0)*(k_1_v_x[j]+2*k_2_v_x[j]+2*k_3_v_x[j]+k_4_v_x[j]);
            v_y[j] = v_y[j] + h*(1.0/6.0)*(k_1_v_y[j]+2*k_2_v_y[j]+2*k_3_v_y[j]+k_4_v_y[j]);
            v_z[j] = v_z[j] + h*(1.0/6.0)*(k_1_v_z[j]+2*k_2_v_z[j]+2*k_3_v_z[j]+k_4_v_z[j]);
        }
        
        for(k=0;k<n_points;k++){
            fprintf(in," %f %f %f %f ", x[k], y[k], v_x[k], v_y[k]);
        }
        
        kinetic = get_kinetic(v_x, v_y, v_z, mass, n_points);
        potential = get_potential(x, y, z, mass, n_points);
        
        fprintf(in,"%f %f \n", kinetic,potential);
    }
    fclose(in);
    
}
     
     */
    
}


void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points, int *ID){
    int i,j;
    FLOAT r_ij;
    for(i=0;i<n_points;i++){
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;
        
        
        /*Las partículas sólo sienten la masa de aquellas que tienen ID negativo*/
        
        for(j=0;j<n_points;j++){
            if(j!=i && ID[j]<0){
                r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
                r_ij = sqrt(r_ij);
                ax[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
                ay[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (y[i] - y[j]);
                az[i] += -G_GRAV *mass[j] / pow(r_ij,1.5) * (z[i] - z[j]);
                
            }
        }
        }
    }
}

FLOAT get_kinetic(FLOAT *v_x, FLOAT *v_y, FLOAT *v_z, FLOAT *mass, int n_points){
    int i;
    FLOAT kinetic;
    kinetic = 0.0;
    
    
    for (i=0; i<n_points; i++) {
        kinetic += 0.5*mass[i]*(pow(v_x[i],2)+pow(v_y[i],2)+pow(v_z[i],2));
    }
    return kinetic;
}

FLOAT get_potential(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points){
    int i,j;
    FLOAT r_ij;
    FLOAT potential;
    potential = 0.0;
    for (i=0; i<n_points; i++) {
        for (j=0; j<n_points; j++) {
            if (i!=j) {
                r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
                r_ij = sqrt(r_ij);
                potential += -G_GRAV * (mass[i] * mass[j]) / r_ij;
            }
        }
    }
    return potential/2;
}

FLOAT * get_memory(int n_points){
    FLOAT * x; 
    if(!(x = malloc(sizeof(FLOAT) * n_points))){
        printf("problem with memory allocation");
        exit(1);
    }
    return x;
}