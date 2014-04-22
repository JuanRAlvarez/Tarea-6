import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Importar datos
data = np.loadtxt("out.txt")


#Separar x e y
ind = [2*i for i in xrange(0 , np.shape(data)[1]/2)]
x = data[:,ind]
ind = [2*i + 1 for i in xrange(0 , np.shape(data)[1]/2)]
y = data[:,ind]

#Preparar graficaas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in xrange(0, np.shape(x)[0]):
    z = (np.shape(x)[0]-i)*np.ones(np.shape(x)[1])
    ax.plot(x[i,:], y[i,:], z)
ax.set_xlabel("Numero de presas")
ax.set_ylabel("Numero de casadores")
ax.set_zlabel("Numero inicial de presas")
ax.set_title("Diferentes curvas solucion parametricas")

mpl.pyplot.savefig("plot_solutions.pdf")
