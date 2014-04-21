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
  FLOAT r = 100.0;//initial orbit radius
  int i;
  for( i=1 ; i<=n_points ; i++){
    z[i] = 0;
    x[i] = r*cos(i*TAU/3);
    y[i] = r*sin(i*TAU/3);
  }
}



/*
 *Initializes velocities
 */
void initvel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points){
  FLOAT r = 100.0;// initial orbit radius
  FLOAT v_0 = sqrt(11*G_GRAV/(3*r));
  int i;
  for( i=1 ; i<=n_points ; i++){
    vz[i] = 0;
    vx[i] = - v_0 * sin(i*TAU/3);
    vy[i] = v_0 * cos(i*TAU/3);
  }
}



/*
 *Initializes accelerations
 */
void initacc(FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points){
  
  FLOAT r = 100.0; //initial orbit radius
  /* Initial accelaration v**2/r */
  FLOAT v_0sq = 11*G_GRAV/(3*r*r);
  int i;

  for( i=1 ; i<=n_points ; i++ ){
    az[i] = 0;
    ax[i] = - v_0sq * cos(i*TAU/3);
    ay[i] = - v_0sq * sin(i*TAU/3);
  } 
}



/*
 *Initializes masses
 */
void initmass(FLOAT *mass, int n_points){
  int i;
  for( i=1 ; i<=n_points ; i++){
    mass[i] = i;
  }
}

/*
 *Calculates the force on the ith partice due to the others
 */

float *force(int i, FLOAT *x , FLOAT *y, float *z , FLOAT * mass, int n_points){

  int j;
  FLOAT F_x = 0;
  FLOAT F_y = 0;
  FLOAT F_z = 0;
  FLOAT *F = malloc(3*sizeof(FLOAT));

  for( j = 1 ; j<=n_points ; j++){
    if( i != j ){
      float epsilon = 0.0; //modificacion a la gravedad
      FLOAT norm_cube = 
	pow( pow(x[i]-x[j],2.0) + pow(y[i]-y[j],2.0) + 
	     pow(z[i]-z[j],2.0) +epsilon, 3.0/2.0);
      F_x = F_x - G_GRAV*mass[i]*mass[j]*(x[i] - x[j])/norm_cube;
      F_y = F_y - G_GRAV*mass[i]*mass[j]*(y[i] - y[j])/norm_cube;
      F_z = F_z - G_GRAV*mass[i]*mass[j]*(z[i] - z[j])/norm_cube;
    }
  } 

  F[0] = F_x;
  F[1] = F_y;
  F[2] = F_z;

  return F;

}

/*
 *Performs RK4
 */
void RK4step(FLOAT *x0 , FLOAT *y0 , FLOAT *z0 , 
	FLOAT *vx0 , FLOAT *vy0 , FLOAT *vz0 , 
	FLOAT *mass, FLOAT delta_t , int n_points){

  int i;

  FLOAT *k1_x, *k1_y, *k1_z;
  FLOAT *k2_x, *k2_y, *k2_z;
  FLOAT *k3_x, *k3_y, *k3_z;
  FLOAT *k4_x, *k4_y, *k4_z;

  k1_x = getmemory(n_points);
  k1_y = getmemory(n_points);
  k1_z = getmemory(n_points);

  k2_x = getmemory(n_points);
  k2_y = getmemory(n_points);
  k2_z = getmemory(n_points);

  k3_x = getmemory(n_points);
  k3_y = getmemory(n_points);
  k3_z = getmemory(n_points);

  k4_x = getmemory(n_points);
  k4_y = getmemory(n_points);
  k4_z = getmemory(n_points);

  FLOAT *k1_vx, *k1_vy, *k1_vz;
  FLOAT *k2_vx, *k2_vy, *k2_vz;
  FLOAT *k3_vx, *k3_vy, *k3_vz;
  FLOAT *k4_vx, *k4_vy, *k4_vz;

  k1_vx = getmemory(n_points);
  k1_vy = getmemory(n_points);
  k1_vz = getmemory(n_points);

  k2_vx = getmemory(n_points);
  k2_vy = getmemory(n_points);
  k2_vz = getmemory(n_points);

  k3_vx = getmemory(n_points);
  k3_vy = getmemory(n_points);
  k3_vz = getmemory(n_points);

  k4_vx = getmemory(n_points);
  k4_vy = getmemory(n_points);
  k4_vz = getmemory(n_points);

  FLOAT *mean_kx, *mean_ky, *mean_kz;
  FLOAT *mean_kvx, *mean_kvy, *mean_kvz;

  mean_kx = getmemory(n_points);
  mean_ky = getmemory(n_points);
  mean_kz = getmemory(n_points);

  mean_kvx = getmemory(n_points);
  mean_kvy = getmemory(n_points);
  mean_kvz = getmemory(n_points);

  
for(i = 1 ; i<=n_points ; i++){
    k1_x[i] = vx0[i];
    k1_y[i] = vy0[i];
    k1_z[i] = vz0[i];

    FLOAT *f = force(i, x0, y0, z0, mass, n_points);
    k1_vx[i] = f[0];
    k1_vy[i] = f[1];
    k1_vz[i] = f[2];
  }
    //first step
 FLOAT *x1, *y1, *z1;
 x1 = getmemory(n_points);
 y1 = getmemory(n_points);
 z1 = getmemory(n_points);

 FLOAT *vx1, *vy1, *vz1;
 vx1 = getmemory(n_points);
 vy1 = getmemory(n_points);
 vz1 = getmemory(n_points);

 for( i=1 ; i<=n_points ; i++){
    x1[i] = x0[i] + (delta_t/2.0) * k1_x[i];
    y1[i] = y0[i] + (delta_t/2.0) * k1_y[i];
    z1[i] = z0[i] + (delta_t/2.0) * k1_z[i];

    vx1[i] = vx0[i] + (delta_t/2.0) * k1_vx[i];
    vy1[i] = vy0[i] + (delta_t/2.0) * k1_vy[i];
    vz1[i] = vz0[i] + (delta_t/2.0) * k1_vz[i];

    k2_x[i] = vx1[i];
    k2_y[i] = vy1[i];
    k2_z[i] = vz1[i];

    FLOAT *f = force(i,x1,y1,z1,mass,n_points);
    k2_vx[i] = f[0]/mass[i];
    k2_vy[i] = f[1]/mass[i];
    k2_vz[i] = f[2]/mass[i];
 }

 //second step
 FLOAT *x2, *y2, *z2;
 x2 = getmemory(n_points);
 y2 = getmemory(n_points);
 z2 = getmemory(n_points);

 FLOAT *vx2, *vy2, *vz2;
 vx2 = getmemory(n_points);
 vy2 = getmemory(n_points);
 vz2 = getmemory(n_points);

 for(i = 1; i<=n_points ; i++){
    x2[i] = x0[i] + (delta_t/2.0) * k2_x[i];
    y2[i] = y0[i] + (delta_t/2.0) * k2_y[i];
    z2[i] = z0[i] + (delta_t/2.0) * k2_z[i];

    vx2[i] = vx0[i] + (delta_t/2.0) * k2_vx[i];
    vy2[i] = vy0[i] + (delta_t/2.0) * k2_vy[i];
    vz2[i] = vz0[i] + (delta_t/2.0) * k2_vz[i];

    k3_x[i] = vx2[i];
    k3_y[i] = vy2[i];
    k3_z[i] = vz2[i];

    FLOAT *f = force(i,x2,y2,z2,mass,n_points);
    k3_vx[i] = f[0]/mass[i];
    k3_vy[i] = f[1]/mass[i];
    k3_vz[i] = f[2]/mass[i];
 }
 
   //third step
    
 FLOAT *x3, *y3, *z3;
 x3 = getmemory(n_points);
 y3 = getmemory(n_points);
 z3 = getmemory(n_points);

 FLOAT *vx3, *vy3, *vz3;
 vx3 = getmemory(n_points);
 vy3 = getmemory(n_points);
 vz3 = getmemory(n_points);

 for( i=1 ; i<=n_points ; i++){ 
    x3[i] = x0[i] + delta_t * k3_x[i];
    y3[i] = y0[i] + delta_t * k3_y[i];
    z3[i] = z0[i] + delta_t * k3_z[i];

    vx3[i]= vx0[i] + delta_t * k2_vx[i];
    vy3[i] = vy0[i] + delta_t * k2_vy[i];
    vz3[i] = vz0[i] + delta_t * k2_vz[i];

    k4_x[i] = vx3[i];
    k4_y[i] = vy3[i];
    k4_z[i] = vz3[i];

    FLOAT *f = force(i, x3, y3, z3,mass, n_points);
    k4_vx[i] = f[0];
    k4_vy[i] = f[1];
    k4_vz[i] = f[2];
 }


    //fouth step
 for( i=0 ; i<=n_points ; i++){
    FLOAT mean_kx = (1/6.0)*(k1_x[i] + 2.0*k2_x[i] + 2.0*k3_x[i] + k4_x[i]);
    FLOAT mean_ky = (1/6.0)*(k1_y[i] + 2.0*k2_y[i] + 2.0*k3_y[i] + k4_y[i]);
    FLOAT mean_kz = (1/6.0)*(k1_z[i] + 2.0*k2_z[i] + 2.0*k3_z[i] + k4_z[i]);

    FLOAT mean_kvx = (1/6.0)*(k1_vx[i] + 2.0*k2_vx[i] + 2.0*k3_vx[i] + k4_vx[i]);
    FLOAT mean_kvy = (1/6.0)*(k1_vy[i] + 2.0*k2_vy[i] + 2.0*k3_vy[i] + k4_vy[i]);
    FLOAT mean_kvz = (1/6.0)*(k1_vz[i] + 2.0*k2_vz[i] + 2.0*k3_vz[i] + k4_vz[i]);

    //New values
    FLOAT x_new = x0[i] + delta_t*mean_kx;
    FLOAT y_new = y0[i] + delta_t*mean_ky;
    FLOAT z_new = z0[i] + delta_t*mean_kz;

    FLOAT vx_new = vx0[i] + delta_t*mean_kvx;
    FLOAT vy_new = vy0[i] + delta_t*mean_kvy;
    FLOAT vz_new = vz0[i] + delta_t*mean_kvz;

    x0[i] = x_new;
    y0[i] = y_new;
    z0[i] = z_new;

    vx0[i] = vx_new;
    vy0[i] = vy_new;
    vz0[i] = vz_new;
 }

}

/*
 * Simulation starts from t=0 
 */
int main(){

  //constants
  int n_points = 3;
  FLOAT delta_t = 0.2;
  FLOAT max_t = 5000;

  //positions
  FLOAT *x = getmemory(n_points);
  FLOAT *y = getmemory(n_points);
  FLOAT *z = getmemory(n_points);  

  //velocities
  FLOAT *vx = getmemory(n_points);
  FLOAT *vy = getmemory(n_points);
  FLOAT *vz = getmemory(n_points);  

  //accelerations
  FLOAT *ax = getmemory(n_points);
  FLOAT *ay = getmemory(n_points);
  FLOAT *az = getmemory(n_points);  

  //mass
  FLOAT *mass = getmemory(n_points);

  //Initialization
  initpos( x , y , z , n_points );
  initvel( vx , vy , vz , n_points );
  initacc( ax , ay , az , n_points );
  initmass( mass , n_points );
  

  //Perform simulation using RK4
  FLOAT t = 0;
  int i;
  printf("%c %c %c %s %s %s\n", 'x', 'y', 'z',"vx", "vy", "vz");
  int j = 0;
  int n_f = 100;
  while( t <= max_t ){

    if(j%n_f == 0){
      for( i=1 ; i <= n_points ; i++){
	printf("%f %f %f %f %f %f\n", x[i] , y[i], z[i], vx[i] , vy[i], vz[i]);
      }
    }

      RK4step( x , y , z ,
	       vx , vy , vz ,
	       ax , ay , az ,
	       mass, delta_t ,n_points);
      t = t + delta_t;
    }
  
  return 0;
}
