all : plot_nbody.dat nbody.x plot_nbody.pdf

plot_nbody.dat : nbody.x
	./nbody.x > plot_nbody.dat
nbody.x : nbody.c
	cc nbody.c -o nbody.x -lm
plot_nbody.pdf : plot_nbody.py
	python plot_nbody.py
clean:
	rm nbody.x plot_nbody.pdf plot_nbody.dat
