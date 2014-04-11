#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TAU (6.28318531)
#define FLOAT float
#define G_GRAV (39.486) 

/*
 * Units:
 * Mass: solar masses
 * Distance: astronomical units
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

}



/*
 *Initializes velocities
 */
void initvel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points){

}



/*
 *Initializes masses
 */
void initmass(FLOAT *mass, int n_points){

}

/*
 * Calculates the force on the ith partice due to the others using
 * fast multipole methods
 */

float *force(int i, FLOAT *x , FLOAT *y, float *z , FLOAT * mass, int n_points){


}


/*
 * Performs integration step using ???
 */
float  integrateStep(FLOAT *x0 , FLOAT *y0 , FLOAT *z0 , 
	FLOAT *vx0 , FLOAT *vy0 , FLOAT *vz0 ,  
	FLOAT *mass, int n_points){
 

}

/*
 * Where the magic happens*.*
 */
int main(){

  //constants
  int n_points;
  FLOAT max_t;

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
    t = integrateStep(x , y , z , x , vy , vz , 
		      mass, n_points);
    //printf("%c %c %c %s %s %s\n", 'x', 'y', 'z',"vx", "vy", "vz");
    // integration method undecided
    }
  
  return 0;
}
