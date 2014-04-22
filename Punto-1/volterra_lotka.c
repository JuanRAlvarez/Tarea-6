#include <stdio.h>
#include <stdlib.h>
#define FLOAT float

/*Pro tip: es mas facil leer este ccodigo de abajo para arriba*/

/*-----------------------------------------------
 *Punto 1, Lotka-Volterra
 *-----------------------------------------------
 */

/******************************************
 * Apartar memoria
 ******************************************/
FLOAT *getmemory(int n_points){
  return malloc( n_points*sizeof(FLOAT) );
}


/******************************************
 * Derivada de x
 ******************************************/

FLOAT dx(FLOAT x , FLOAT y){
  FLOAT a , b;
  a  =20.0;
  b = 1.0;
  return a*x - b*x*y ; 
}


/******************************************
 * Derivada de y
 ******************************************/

FLOAT dy(FLOAT x , FLOAT y){
  FLOAT c , d;
  c  =30.0;
  d = 1.0;
  return -c*y + d*x*y ; 
}


/*******************************************
 *-------------Runge-Kutta 4 step-----------
 *******************************************/
FLOAT *RK4step(FLOAT x0 , FLOAT y0 , FLOAT delta_t){

  FLOAT *xy = getmemory(2);
  FLOAT x, y;

  FLOAT k1x, k2x, k3x, k4x;
  FLOAT k1y, k2y, k3y, k4y;

  k1x = delta_t*dx(x0,y0);
  k1y = delta_t*dy(x0,y0);

  k2x = delta_t*dx(x0 + k1x/2.0 , y0 + k1y/2.0);
  k2y = delta_t*dy(x0 + k1x/2.0 , y0 + k1y/2.0);

  k3x = delta_t*dx(x0 + k2x/2.0 , y0 + k2y/2.0);
  k3y = delta_t*dy(x0 + k2x/2.0 , y0 + k2y/2.0);

  k4x = delta_t*dx(x0 + k3x , y0 + k3y);
  k4y = delta_t*dy(x0 + k3x , y0 + k3y);

  x = x0 + (1/6.0)*(k1x + k2x/2.0 + k3x/2.0 + k4x);
  y = y0 + (1/6.0)*(k1y + +k2y/2.0 + k3y/2.0 + k4y);

  /* En este caso la funcion es siempre mayor o igual a 0 */
  if(x >= 0){
    xy[0] = x;
  }
  else{
    xy[0] = 0;
  }

  if(y >= 0 ){
    xy[1] = y;
  }
  else{
    xy[1] = 0;
  }
  return xy;

}

/*******************************************
 *-------------Runge-Kutta 4----------------
 *******************************************/

FLOAT **RK4(FLOAT x0 , FLOAT y0 , FLOAT init_t , 
	    FLOAT delta_t , int n_points ){
  FLOAT **xy;
  xy = malloc(2*sizeof(FLOAT*));
  FLOAT *x = getmemory(n_points);
  FLOAT *y = getmemory(n_points);
  x[0] = x0;
  y[0] = y0;
  FLOAT *xyi = getmemory(2);
  int i;
  for( i = 1 ; i<=n_points ; i++ ){
    xyi = RK4step( x[i-1] , y[i-1] , delta_t );
    x[i] = xyi[0];
    y[i] = xyi[1];
  }
  xy[0] = x;
  xy[1] = y;
  return xy;
  }


/************************************************
 *Resolver el punto 1
 *Escribe los datos en un archivo de texto de la 
 *siguiente forma:
 *cada linea es una serie de tiempo de la forma 
 *x0 y0 x1 y1 x2 y2 x2 y3 ...
 ************************************************/
void punto1(FLOAT x0 , int x0i , FLOAT y0,  FLOAT init_t , FLOAT delta_t , int n_points){

  /* Disminuir x0 uno a uno hasta uno*/

  FILE *f;
  f = fopen("out.txt" , "w");
  int i;
  for( ; x0i >= 1 ; x0i-- ){
    FLOAT **xy = RK4( x0 , y0 , init_t , delta_t , n_points );
    FLOAT *x = xy[0];
    FLOAT *y = xy[1];
    for( i = 0 ; i < n_points ; i++){
      fprintf( f , "%f %f " , x[i] , y[i]);
    } 
    fprintf( f , "\n");
    x0 = x0 - 1;
  }
  fclose(f);

}

/***************************************************
 *El valor de equilirio para x es C/D  = 30
 *El valor de equilibrio para y es A/B = 20
 ***************************************************/

int main(){

  /*
   *Variables
   */
  FLOAT delta_t = 0.0001; //Potencia de 10
  FLOAT init_t = 0;
  FLOAT final_t = 1;
  int  n_points = (int)( ( final_t - init_t ) / delta_t );
  FLOAT x0 = 30.0;
  FLOAT y0 = 20.0;
  int x0i = (int)x0;

/* Resolver el punto  */
  punto1(x0 , x0i , y0,  init_t , delta_t, n_points);

/* \o/ ¡festejar! \o/  */
  return 0;
}
