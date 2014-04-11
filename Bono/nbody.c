#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TAU (6.28318531)
#define FLOAT float
#define G_GRAV (39.486) 

/*
 * Units:
 * Mass: solar masses
 * Distance: parsecs
 * Time: years  
 */


/*
 *allocates memory
 */
FLOAT *getmemory(int n_points){
  return malloc( n_points*sizeof(FLOAT) );
}


/*
 *Initializes position
 */
void initpos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points){

  int i; 
  FLOAT R = 20;
  FLOAT theta, phi, r, xi, yi, zi;
  for( i = 1; i <= n_points ; i++){
    theta = acos(drand48()*2-1);
    phi = drand48()*TAU;
    r = pow(drand48() , 1/3.0)*R;
    xi = r*sin(theta)*cos(phi);
    yi = r*sin(theta)*sin(phi);
    zi = r*cos(theta);
    x[i] = xi;
    y[i] = yi;
    z[i] = zi; 

  }
}



/*
 *Initializes velocities
 */
void initvel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points){

  int i;
  for(i = 1 ; i<= n_points ; i++){
    vx[i]=0;
    vy[i]=0;
    vz[i]=0;
  }

}



/*
 *Initializes masses
 */
void initmass(FLOAT *mass, int n_points){
    int i;
  for(i = 1 ; i<= n_points ; i++){
    mass[i] = 1;
  }
}

/*
 * Calculates the force per unit acceleration on the ith 
 * partice due to the others using
 * fast multipole methods
 * RETURN: a 3-position pointer acc with
 * ax = acc[0]
 * ay = acc[1]
 * az = acc[2]
 */

float *accel(int i, FLOAT *x , FLOAT *y, float *z , FLOAT *mass, int n_points){



}


/*
 * Performs integration step using the Stormer-Verlet
 * (leapfrog)  method
 */
void  integrateStep(FLOAT *x0 , FLOAT *y0 , FLOAT *z0 , 
		     FLOAT *vx0 , FLOAT *vy0 , FLOAT *vz0 ,
		     FLOAT *mass, FLOAT delta_t, int n_points){

  FLOAT *x, *y, *z, *vx, *vy, *vz; 
  //positions
  x = getmemory(n_points);
  y = getmemory(n_points);
  z = getmemory(n_points);  

  //not necessary to store velocities

  //accelerations
  FLOAT *ax0, *ay0, *az0;
  ax0 = getmemory(n_points);
  ay0 = getmemory(n_points);
  az0 = getmemory(n_points);
  
  FLOAT *acc;
  acc = getmemory(3);

  //La Magia
  int i;
  for( i = 0 ; i<=n_points ; i++ ){
    acc = accel(i, x0, y0, z0, mass, n_points);
    ax0[i] = acc[0];
    ay0[i] = acc[1];
    az0[i] = acc[2];
    x[i] = x0[i] + delta_t*( vx0[i] + acc[0]*(delta_t/2.0) );
    y[i] = y0[i] + delta_t*( vy0[i] + acc[1]*(delta_t/2.0) );
    z[i] = z0[i] + delta_t*( vz0[i] + acc[1]*(delta_t/2.0) );
  } 
  for( i = 0 ; i<=n_points ; i++ ){
    acc = accel(i, x, y, z, mass, n_points);
    vx0[i] = vx0[i] + (delta_t/2.0)*(ax0[i] + acc[0]);
    vy0[i] = vy0[i] + (delta_t/2.0)*(ay0[i] + acc[1]);
    vx0[i] = vz0[i] + (delta_t/2.0)*(az0[i] + acc[2]);
  }
}

/*
 * Where the magic happens*.*
 */
int main(){

  //constants
  int n_points;
  FLOAT max_t;
  FLOAT delta_t;

  //positions
  FLOAT *x = getmemory(n_points);
  FLOAT *y = getmemory(n_points);
  FLOAT *z = getmemory(n_points);  

  //velocities
  FLOAT *vx = getmemory(n_points);
  FLOAT *vy = getmemory(n_points);
  FLOAT *vz = getmemory(n_points);  

  //mass
  FLOAT *mass = getmemory(n_points);

  //Initialization
  initpos( x , y , z , n_points );
  initvel( vx , vy , vz , n_points );
  initmass( mass , n_points );
  

  //Perform simulation using RK4
  FLOAT t = 0;
  int i;

  while( t <= max_t ){
    integrateStep(x , y , z , x , vy , vz , delta_t, 
		  mass, n_points);
    t = t+delta_t;
    //printf("%c %c %c %s %s %s\n", 'x', 'y', 'z',"vx", "vy", "vz");
    // integration method undecided
    }
  
  return 0;
}
