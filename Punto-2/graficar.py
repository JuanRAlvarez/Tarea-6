import pylab as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt("3cuerpos.dat")

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(data[:,0],data[:,1],data[:,2])
ax.plot(data[:,3],data[:,4],data[:,5])
ax.plot(data[:,6],data[:,7],data[:,8])
plt.show()
plt.savefig("Grafica3d.pdf")
