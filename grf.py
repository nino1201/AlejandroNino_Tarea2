import numpy as np
import matplotlib.pyplot as plt

datos=np.genfromtxt('data.dat')
x=datos[:,0]
rho=datos[:,1]
p=datos[:,2]
v=datos[:,3]
plt.subplot(311)
plt.scatter(x,p)
plt.title('Presion')

plt.subplot(312)
plt.scatter(x,rho)
plt.title('Densidad')

plt.subplot(313)
plt.scatter(x,v)
plt.title('Velocidad')
plt.show()
