import numpy as np
import matplotlib as plt

datos = np.loadtxt("IC.dat")

N = np.shape(datos)[0]

for i in range(N):
    plt.scatter(datos[i,0],datos[i,1])


plt.show()