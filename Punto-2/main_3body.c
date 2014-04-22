#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1

FLOAT * get_memory(int n_points);
void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius);
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius);
void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);
FLOAT get_kinetic(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *mass, int n_points);
FLOAT get_potential(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *mass, int n_points);

/* ESTE CÓDIGO ESTÁ BASADO EN EL CÓDIGO DEL PROFESOR JAIME FORERO (github.com/forero), PERO TIENE ALGUNAS MODIFICACIONES, A SABER:
 se implementa el método de runge kutta de cuarto orden
 */

int main(int argc, char **argv){
  
  /*positions of all particles*/
  FLOAT *x;
  FLOAT *y;
  FLOAT *z;
  
  /*velocities of all particles*/
  FLOAT *v_x;
  FLOAT *v_y;
  FLOAT *v_z;

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
    
    
  /*masses*/
  FLOAT *mass;

  /*timestep variables*/
  FLOAT h= 0.001;
  int n_steps = (int)(100/h);
  int n_points = 3;
  FLOAT radius = 100.0;
  FLOAT unit_mass = 1.0; 
  FLOAT vel_initial = sqrt((11.0/3.0) * G_GRAV * unit_mass / (sqrt(3.0)*radius));
  FLOAT kinetic;
  FLOAT potential;
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

  initialize_pos(x,y,z, n_points, radius);
  initialize_vel(v_x,v_y,v_z, n_points, vel_initial, radius);
  initialize_mass(mass, n_points, unit_mass);

  /*implementation of a second order runge kutta integration*/
    
    FILE *in;
    in = fopen("3cuerpos.dat","w");
    
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
          fprintf(in," %f %f %f ", x[k], y[k], z[k]);
      }
      kinetic = get_kinetic(x, y, z, v_x, v_y, v_z, a_x, a_y, a_z, mass, n_points);
      potential = get_potential(x, y, z, v_x, v_y, v_z, a_x, a_y, a_z, mass, n_points);
      fprintf(in,"%f %f \n",kinetic, potential);
  }
    fclose(in);
    
}

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points){
  int i,j;
  FLOAT r_ij;
  for(i=0;i<n_points;i++){
    ax[i]=0.0;
    ay[i]=0.0;
    az[i]=0.0;
    
    for(j=0;j<n_points;j++){
      if(j!=i){
          r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
          r_ij = sqrt(r_ij);
          ax[i] += -G_GRAV *mass[j]/ pow(r_ij,3) * (x[i] - x[j]);
          ay[i] += -G_GRAV *mass[j]/ pow(r_ij,3) * (y[i] - y[j]);
          az[i] += -G_GRAV *mass[j] / pow(r_ij,3) * (z[i] - z[j]);
          
      }
    }    
  }  
}

void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;
  
  for(i=0;i<n_points;i++){
    x[i] = cos(delta_theta * i) * radius;
    y[i] = sin(delta_theta * i) * radius;
    z[i] = 0.0;
  }
}

void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;
  
  for(i=0;i<3;i++){
   // vx[i] = -sin(delta_theta * i) * vel;
   // vy[i] = cos(delta_theta * i) * vel;
      vx[i] = 0.0;
      vy[i] = 0.0;
      vz[i] = vel;
  }  

    
}

void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass){
  int i;
  for (i=0;i<n_points;i++){
    mass[i] = (i+1) * unit_mass;
  }
}

FLOAT get_kinetic(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *mass, int n_points){
    
    FLOAT kinetic;
    int i;
    for (i=0; i<n_points; i++) {
        kinetic += mass[i]*(pow(vx[i],2)+pow(vy[i],2)+pow(vz[i],2));
    }
    
    return kinetic;
}

FLOAT get_potential(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *mass, int n_points){
    
    FLOAT potential;
    FLOAT r_ij;
    int i,j;
    for (i=0; i<n_points; i++) {
        for (j=i+1; j<n_points; j++) {
                r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
                r_ij = sqrt(r_ij);
                potential += G_GRAV*mass[i]*mass[j]/(r_ij);
            
        }
    }
    
    return -potential;
}

FLOAT * get_memory(int n_points){
  FLOAT * x; 
  if(!(x = malloc(sizeof(FLOAT) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}


