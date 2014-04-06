import pylab as plt
import numpy as np

datos = np.loadtxt("3cuerpos.dat")

plt.scatter(datos[:,0],datos[:,1],c='r')
plt.scatter(datos[:,2],datos[:,3],c='g')
plt.scatter(datos[:,4],datos[:,5],c='b')

plt.show()