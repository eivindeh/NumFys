import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
file = open('EulerOut.txt')

data = np.loadtxt('EulerOut.txt')

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(data[:,0],data[:,1],data[:,2])
plt.show()
