all: Grafica.pdf

Grafica.pdf: 1.dat 2.dat 3.dat 4.dat 5.dat plots.py
	python plots.py 5

1.dat 2.dat 3.dat 4.dat 5.dat: evoluciona.x Conds.dat
	./evoluciona.x 5 5 Conds.dat

evoluciona.x: evoluciona.c Conds.dat
	cc evoluciona.c -lm -o evoluciona.c

Conds.dat: IC1.dat IC2.dat IC.x
	cat IC1.dat IC2.dat >> Conds.dat

IC1.dat: IC.x
	./IC.x 0 0 0 0 0 0 24 1 20 1000 IC1.dat

IC2.dat: IC.x
	./IC.x 150 200 0 0 -100 0 24 1 20 100 IC2.dat

IC.x: IC.c
	cc IC.c -lm -o IC.x

clean:
	rm IC.x
	rm IC1.dat
	rm IC2.dat
	rm 1.dat
	rm 2.dat
	rm 3.dat
	rm 4.dat
	rm 5.dat
	rm evoluciona.x
	
	