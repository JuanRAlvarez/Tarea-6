#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

    float delta_t;
    int N;
    
    reservar.memoria();
    inicializar_pos();
    inicializar_vel();
    inicializar_masas();
    
    for (i=0; i<N; i++) {
        update_pos();
        update_vel();

    }

}