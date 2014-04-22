import numpy as np
import matplotlib.pyplot as plt
import sys

l = int(sys.argv[1])

for i in range(l):
    nombre = str(i+1)+'.dat'
    print nombre
    datos = np.loadtxt(nombre)
    N = np.shape(datos)[0]

    for i in range(N):
        if (datos[i,0] == -1):
           plt.scatter(datos[i,1],datos[i,2],c='r')
        else:
           plt.scatter(datos[i,1],datos[i,2],c='g')
    plt.show()