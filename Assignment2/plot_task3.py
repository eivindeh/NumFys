import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
file = open('EulerOut.txt')

data = np.loadtxt('EulerOut.txt')

fig = plt.figure()
ax = fig.gca(projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b')

ax.set_xlim3d(-4,4)
ax.set_ylim3d(-4,4)
ax.set_zlim3d(-4,4)

ax.plot(data[:,0],data[:,1],data[:,2])
plt.show()
