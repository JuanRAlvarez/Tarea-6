/*terms for runge kutta*/
FLOAT *xmed;
FLOAT *ymed;
FLOAT *zmed;
FLOAT *xold;
FLOAT *yold;
FLOAT *zold;
FLOAT *x1;
FLOAT *x2;
FLOAT *x3;
FLOAT *y1;
FLOAT *y2;
FLOAT *y3;
FLOAT *v_x1;
FLOAT *v_x2;
FLOAT *v_x3;
FLOAT *v_y1;
FLOAT *v_y2;
FLOAT *v_y3;
FLOAT *v_z1;
FLOAT *v_z2;
FLOAT *v_z3;
FLOAT *v_xmed;
FLOAT *v_ymed;
FLOAT *v_zmed;
FLOAT *v_xold;
FLOAT *v_yold;
FLOAT *v_zold;
FLOAT *a_xold;
FLOAT *a_yold;
FLOAT *a_zold;
FLOAT *k1_x;
FLOAT *k1_p_x;
FLOAT *k2_x;
FLOAT *k2_p_x;
FLOAT *k3_x;
FLOAT *k3_p_x;
FLOAT *k4_x;
FLOAT *k4_p_x;
FLOAT *k1_y;
FLOAT *k1_p_y;
FLOAT *k2_y;
FLOAT *k2_p_y;
FLOAT *k3_y;
FLOAT *k3_p_y;
FLOAT *k4_y;
FLOAT *k4_p_y;
FLOAT *k1_z;
FLOAT *k1_p_z;
FLOAT *k2_z;
FLOAT *k2_p_z;
FLOAT *k3_z;
FLOAT *k3_p_z;
FLOAT *k4_z;
FLOAT *k4_p_z;


/*memory allocation*/


xold[i] = x[i];
v_xold[i] = v_x[i];
a_xold[i] = a_x[i];
yold[i] = y[i];
v_yold[i] = v_y[i];
a_yold[i] = a_y[i];
zold[i] = z[i];
v_zold[i] = v_z[i];
a_zold[i] = a_z[i];


//VALORES DE K
k1_x[i] = xold[i];
k1_y[i] = yold[i];
k1_z[i] = zold[i];

k1_p_x[i] = v_xold[i];
k1_p_y[i] = v_yold[i];
k1_p_z[i] = v_zold[i];


//primer paso

x1[i] = xold[i] + (h/2.0)*k1_x[i];
y1[i] = yold[i] + (h/2.0)*k1_y[i];
z1[i] = zold[i] + (h/2.0)*k1_z[i];

v_x1[i] = v_xold[i] + (h/2.0)*k1_p_x[i];
v_y1[i] = v_yold[i] + (h/2.0)*k1_p_y[i];
v_z1[i] = v_zold[i] + (h/2.0)*k1_p_z[i];

k2_x[i]= v_x1[i];
k2_p_x[i] = a_x


/* k_1_prime1 = func_prime_1(x_old,y1_old, y2_old)
 k_1_prime2 = func_prime_2(x_old,y1_old, y2_old)
 
 #first step
 x1 = x_old+ (h/2.0)
 y1_1 = y1_old + (h/2.0) * k_1_prime1
 y2_1 = y2_old + (h/2.0) * k_1_prime2
 k_2_prime1 = func_prime_1(x1, y1_1, y2_1)
 k_2_prime2 = func_prime_2(x1, y1_1, y2_1)
 
 #second step
 x2 = x_old + (h/2.0)
 y1_2 = y1_old + (h/2.0) * k_2_prime1
 y2_2 = y2_old + (h/2.0) * k_2_prime2
 k_3_prime1 = func_prime_1(x2, y1_2, y2_2)
 k_3_prime2 = func_prime_2(x2, y1_2, y2_2)
 
 
 #third
 x3 = x_old + h
 y1_3 = y1_old + h * k_3_prime1
 y2_3 = y2_old + h * k_3_prime2
 k_4_prime1 = func_prime_1(x3, y1_3, y2_3)
 k_4_prime2 = func_prime_2(x3, y1_3, y2_3)
 
 #fourth step
 average_k_1 = (1.0/6.0)*(k_1_prime1 + 2.0*k_2_prime1 + 2.0*k_3_prime1 + k_4_prime1)
 average_k_2 = (1.0/6.0)*(k_1_prime2 + 2.0*k_2_prime2 + 2.0*k_3_prime2 + k_4_prime2)
 
 x_new = x_old + h
 y_1_new = y1_old + h * average_k_1
 y_2_new= y2_old + h * average_k_2
 return x_new, y_1_new, y_2_new*/

//Middle point
xmed[j] = x[j] + (h/2.0)*v_x[j];
ymed[j] = y[j] + (h/2.0)*v_y[j];
zmed[j] = z[j] + (h/2.0)*v_z[j];