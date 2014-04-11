all: Grafica3d.pdf

Grafica3d.pdf: 3cuerpos.dat graficar.py
	python graficar.py

3cuerpos.dat: 3cuerpos.x
	./3cuerpos.x

3cuerpos.x: main_3body.c
	cc main_3body.c -lm -o 3cuerpos.x

clean:
	rm 3cuerpos.x
	rm 3cuerpos.dat