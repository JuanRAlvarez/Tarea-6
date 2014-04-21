#include<stdio.h>
#include<stdlib.h>

void escribir(float *x, float *y, float *z, int pasos, char *filename);

int main(int argc, char **argv){
    
    float lista1[10];
    float lista2[10];
    float lista3[10];

    int i;
    
    
    for (i=0; i<10; i++) {
        lista1[i] = i*i*i;
        lista2[i] = i*i;
        lista3[i] = i*3.0;
        }
    
    escribir(lista1,lista2,lista3, 10, "hola.dat");
    
    return 0;
}

void escribir(float *x, float *y, float *z, int pasos, char *filename){

    
    FILE *in;
    
    in = fopen(filename,"w");
    int i;
    
    for (i=0; i<pasos; i++) {
        fprintf(in,"%f %f %f \n", x[i], y[i], z[i]);
    }
    
    fclose(in);

}
