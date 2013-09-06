#!/usr/bin/python

import numpy as np
import ConfigParser
import matplotlib.pyplot as plt

def build_area():
    area = ConfigParser.RawConfigParser()
    area.read('area.config')

    nx = int(area.get('Area', 'Nx'))
    ny = int(area.get('Area', 'Ny'))
    nz = int(area.get('Area', 'dNz'))

    if nx != ny:
        raise Exception("This method expects that nx == ny!")

    lx = float(area.get('Area', 'Lx'))
    ly = float(area.get('Area', 'Ly'))
    lz = float(area.get('Area', 'Lz'))

    R = float(area.get('Circle', 'R'))
    z0 = float(area.get('Circle', 'x0'))
    y0 = float(area.get('Circle', 'y0'))

    if R is None or z0 is None or y0 is None:
        raise Exception("There is no circle parameters!")

    x = np.linspace(0.01, lx-0.01, nx)
    y = np.linspace(0.01, ly-0.01, ny)
    z = np.linspace(0.01, lz-0.01, nz)

    y1_array = -np.sqrt(R**2 - (z-z0)**2) + z0
    y2_array = np.sqrt(R**2 - (z-z0)**2) + z0
    y2 = np.concatenate((y1_array, y2_array))
    y2 = y2[::2]
    y2 = np.unique(y2)
    y2.sort()
    y = y2

    #z1_array = -np.sqrt(R**2 - (y-y0)**2) + y0
    #z2_array = np.sqrt(R**2 - (y-y0)**2) + y0
    #z2 = np.concatenate((z1_array, z2_array))
    ##z2 = z2[::2]
    #z2 = np.unique(z2)
    #z2.sort()

    #z = np.concatenate((z, z2))
    #y = np.concatenate((y, y2))
    ##z = z[::2]
    #z = np.unique(z)
    #z.sort()
    ##y = y[::2]
    #y = np.unique(y)
    #y.sort()

    #Z, Y = np.meshgrid(z,y)
    #circle = plt.Circle((0.5, 0.5), radius=0.5)
    #figure = plt.gcf()
    #figure.gca().add_artist(circle)
    #plt.plot(Z, Y, 'ro')
    #plt.show()

    np.savetxt('area.x.coord', x, delimiter='\n', fmt='%2.8f')
    np.savetxt('area.y.coord', y, delimiter='\n', fmt='%2.8f')
    np.savetxt('area.z.coord', z, delimiter='\n', fmt='%2.8f')

if __name__ == '__main__':
    build_area()
