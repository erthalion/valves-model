import numpy as np
import matplotlib.pyplot as plt

lx = 1
ly = 1
nx = 50
ny = 50
R = 0.5
x0 = 0.5
y0 = 0.5

x = np.linspace(0, lx, nx)

y1_array = -np.sqrt(R**2 - (x-x0)**2) + x0
y2_array = np.sqrt(R**2 - (x-x0)**2) + x0

y = np.concatenate((y1_array, y2_array))

#debug
#X, Y = np.meshgrid(x,y)

#circle = plt.Circle((0.5, 0.5), radius=0.5)

#figure = plt.gcf()
#figure.gca().add_artist(circle)

#plt.plot(X, Y, 'ro')
#plt.show()
