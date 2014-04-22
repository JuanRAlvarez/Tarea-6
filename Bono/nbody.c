#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TAU (6.28318531)
#define FLOAT float
#define G_GRAV (4499.93) 
#define EPSILON (0.01)



/******************************************
 *BONO 
 *
 * Units:
 * Mass: solar masses
 * Distance: parsecs
 * Time: 10E9 years  
 ******************************************/

/******************************************
 *Separar memoria 
 ******************************************/
FLOAT *getmemory(int n_masses){
  return malloc( n_masses*sizeof(FLOAT) );
}

/******************************************
 *Inicializar posiciones
 ******************************************/
FLOAT **initpos(int n_masses){
  FLOAT **xyz = malloc(3*sizeof(FLOAT*));
  FLOAT *x = getmemory(n_masses);
  FLOAT *y = getmemory(n_masses);
  FLOAT *z = getmemory(n_masses);
  int i; 
  FLOAT R = 20;
  FLOAT theta, phi, r, xi, yi, zi;
  for( i = 1; i <= n_masses ; i++){
    theta = acos(drand48()*2-1);
    phi = drand48()*TAU;
    r = pow(drand48() , 1/3.0)*R;
    xi = r*sin(theta)*cos(phi);
    yi = r*sin(theta)*sin(phi);
    zi = r*cos(theta);
    x[i-1] = xi;
    y[i-1] = yi;
    z[i-1] = zi; 
  }
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
  return xyz;
}
/******************************************
 *Inicializar velocidades
 ******************************************/
FLOAT **initvel(int n_masses){
  FLOAT **vxyz = malloc(3*sizeof(FLOAT*));
  FLOAT *vx = getmemory(n_masses);
  FLOAT *vy = getmemory(n_masses);
  FLOAT *vz = getmemory(n_masses);
  int i; 
  for( i = 0 ; i<n_masses ; i++){
    vx[i] = 0;
    vy[i] = 0;
    vz[i] = 0;
  }
  vxyz[0] = vx;
  vxyz[1] = vy;
  vxyz[2] = vz;
  return vxyz;
}

/*********************************************************
 * Aceleracion de la i-esima particula
 *********************************************************/
FLOAT *accel(int i, FLOAT *x, FLOAT *y, FLOAT *z, int n_masses){
  int j;
  FLOAT ax, ay, az;
  ax = 0;
  ay = 0;
  az = 0;
  FLOAT *acc;
  acc = getmemory(3);
  FLOAT rijsq;

  for( j = 0 ; j<n_masses ; j++){
    if(i!=j){
      rijsq = pow(x[i] - x[j], 2 ) + pow(y[i] - y[j], 2) 
		 + pow(z[i] - z[j], 2 ) + EPSILON;
      ax = ax + (x[j] - x[i])/pow(rijsq , 3.0/2.0);
      ay = ay + (y[j] - y[i])/pow(rijsq , 3.0/2.0);
      az = az + (z[j] - z[i])/pow(rijsq , 3.0/2.0);

    }
  } 

  acc[0] = G_GRAV*ax;
  acc[1] = G_GRAV*ay;
  acc[2] = G_GRAV*az;
  return acc;

}

/*********************************************************
 * Leap Frog (step) (Stormer -Verlet method)
 *********************************************************/
void  leapfrogStep( FLOAT *x0 , FLOAT *y0 , FLOAT*z0 , 
		      FLOAT *vx0 , FLOAT *vy0 , FLOAT *vz0 , 
		    FLOAT delta_t , int n_masses){

  //accelerations
  FLOAT *ax0, *ay0, *az0;
  ax0 = getmemory(n_masses);
  ay0 = getmemory(n_masses);
  az0 = getmemory(n_masses);
  
  FLOAT *acc;
  acc = getmemory(3);

  //La Magia
  int i;
  for( i = 0 ; i<n_masses ; i++ ){
    acc = accel(i, x0, y0, z0, n_masses);
    ax0[i] = acc[0];
    ay0[i] = acc[1];
    az0[i] = acc[2];
    x0[i] = x0[i] + delta_t*( vx0[i] + acc[0]*(delta_t/2.0) );
    y0[i] = y0[i] + delta_t*( vy0[i] + acc[1]*(delta_t/2.0) );
    z0[i] = z0[i] + delta_t*( vz0[i] + acc[2]*(delta_t/2.0) );
  } 
  for( i = 0 ; i<n_masses ; i++ ){
    acc = accel(i, x0, y0, z0, n_masses);
    vx0[i] = vx0[i] + (delta_t/2.0)*(ax0[i] + acc[0]);
    vy0[i] = vy0[i] + (delta_t/2.0)*(ay0[i] + acc[1]);
    vz0[i] = vz0[i] + (delta_t/2.0)*(az0[i] + acc[2]);
  }

}

/*********************************************************
 *La energia potencial
 *********************************************************/
FLOAT potEn(FLOAT *x , FLOAT *y , FLOAT*z , int n_masses){
  int i;
  int j;
  FLOAT rij;
  FLOAT sum;
  for(i = 0 ; i <n_masses ; i++){
    for(j = 0; j<n_masses ; j++){
      if( i!=j ){
	rij = pow( pow(x[i] - x[j], 2 ) + pow(y[i] - y[j], 2 ) 
		 + pow(z[i] - z[j], 2 ) + EPSILON , 1/2.0 );
	sum = sum - 1/rij;
      }
    }
  }
  return sum*G_GRAV;
}

/**********************************************************
 *Dos veces la energia cinetica
 **********************************************************/
FLOAT kin2En(FLOAT *vx , FLOAT *vy , FLOAT*vz , int n_masses){
  FLOAT sum = 0;
  int i;  
  for(i = 0; i<n_masses; i++){
    sum = sum + pow(vx[i],2) + pow(vy[i], 2) + pow(vz[i],2);
  }
  return sum;
}


/***********************************************************
 *Leapfrog!
 ***********************************************************/
void leapfrog(FLOAT *x0 , FLOAT *y0 , FLOAT *z0 , 
	      FLOAT *vx0 , FLOAT *vy0 , FLOAT *vz0 , 
	      FLOAT delta_t , int n_points , int n_masses){

  //archivo
  int i;
  int j;
  FILE *f;
  FILE *f2;
  f = fopen("plot_nbody.dat", "w");
  f2 = fopen("virEn.dat" , "w");
  for( i = 1 ; i <= n_points ; i++ ){
    leapfrogStep( x0 , y0 , z0 , vx0 , vy0 , vz0 , 
		  delta_t , n_masses);

    //Imprimir 20 momentos
    if(i % n_points/20 == 0 ){
      FLOAT U = potEn(x0 , y0 , z0 , n_masses);
      FLOAT K2 = kin2En(vx0 , vy0 , vz0 , n_masses);
      for(j = 0 ; j < n_masses ; j++){
	fprintf(f , "%f %f %f %f %f %f " , x0[j] , y0[j] , z0[j], 
		vx0[j] , vy0[j] , vz0[j]);
      }
      //Imprimir 2k + u y k+u
      fprintf(f, "\n");
      fprintf(f2, "%f %f\n" , K2 +U, K2/2.0 + U);
    }

  }
  fclose(f);
  fclose(f2);
}


int main(){

  /*
   * Variables
   */

  FLOAT delta_t = 0.001; //Potencia de 10
  FLOAT init_t = 0;
  FLOAT final_t = 1;
  int  n_points = (int)( ( final_t - init_t ) / delta_t );
  int n_masses = 10000;
  FLOAT *x0 , *y0, *z0;
  FLOAT *vx0 , *vy0, *vz0;

  /*
   *Inicializaion
   */  
  FLOAT **xyz = initpos(n_masses);
  x0 = xyz[0];
  y0 = xyz[1];
  z0 = xyz[2]; 

  FLOAT **vxyz = initvel(n_masses);
  vx0 = vxyz[0];
  vy0 = vxyz[1];
  vz0 = vxyz[2]; 
  /*
   *Bono, por favor :)
   */
  leapfrog(x0 , y0 , z0 , vx0 , vy0, vz0 ,
	   delta_t , n_points , n_masses);
  return 0;
}
