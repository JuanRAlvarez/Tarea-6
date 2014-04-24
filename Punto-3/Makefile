all: Grafica.pdf

Grafica.pdf: 1.dat plots.py
	python plots.py 5

1.dat: evolve.x Conds.dat
	./evolve.x 5 5 Conds.dat

evolve.x: evolve.c
	cc evolve.c -lm -o evolve.x

Conds.dat: IC1.dat IC2.dat
	cat IC1.dat IC2.dat >> Conds.dat

IC1.dat: IC.x
	./IC.x 0 0 0 0 0 0 24 1 20 1000 IC1.dat

IC2.dat: IC.x
	./IC.x 150 200 0 0 -100 0 24 1 20 1000 IC2.dat

IC.x: IC.c
	cc IC.c -lm -o IC.x

clean:
	rm *.x *.dat
	
	