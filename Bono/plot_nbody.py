#Plots data from nbody.c
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#Importar datos
data = np.loadtxt("plot_nbody.dat")


#Separar x y z
ind = [6*i for i in xrange(0 , np.shape(data)[1]/6)]
x = data[:,ind]
ind = [6*i + 1 for i in xrange(0 , np.shape(data)[1]/6)]
y = data[:,ind]
ind = [6*i + 2 for i in xrange(0 , np.shape(data)[1]/6)]
z = data[:,ind]

#Preparar graficaas
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[19,1:10], y[19,:], z[19,:])
ax.set_xlabel("x(pc)")
ax.set_ylabel("y(pc)")
ax.set_zlabel("z(pc)")
ax.set_title("Posicion cerca del final del los tiempos")

mpl.pyplot.savefig("plot_nbody.pdf")



#Importar datos
data = np.loadtxt("virEn.dat")

G = data[:,1]
E = data[:,2]

plt.plot(G)
plt.plot(E)
ptt.savefig("virEn.pdf")
