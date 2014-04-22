all: plot_solutions.pdf a.out

plot_solutions.pdf : plot_solutions.py
	python plot_solutions.py
plot_solutions.py : a.out
	cc volterra_lotka.c
a.out : volterra_lotka.c
	cc volterra_lotka.c
clean:
	rm a.out plot_solutions.pdf 
