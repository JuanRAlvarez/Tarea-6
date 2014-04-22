all : plot_nbody.dat virEn.dat plot_nbody.pdf virEn.pdf

plot_nbody.dat : a.out
	./a.out
a.out : nbody.c
	cc nbody.c -lm
virEn.dat : a.out
	./a.out
plot_nbody.dat : a.out
	./a.out
plot_nbody.pdf : plot_nbody.py
	python plot_nbody.py
virEn.pdf : plot_nbody.py
	python plot_nbody.py
clean:
	rm plot_nbody.pdf plot_nbody.dat virEn.dat a.out
