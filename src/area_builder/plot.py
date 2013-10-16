#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt('area.x.coord', delimiter='\n')
y = np.loadtxt('area.y.coord', delimiter='\n')

X, Y = np.meshgrid(x, y)

fig = plt.figure()
circle = plt.Circle((0.5, 0.5), radius=0.5)
figure = plt.gcf()
figure.gca().add_artist(circle)
plt.plot(X, Y, 'ro')
plt.show()
