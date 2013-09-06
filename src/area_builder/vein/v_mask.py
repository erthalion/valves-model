#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import ConfigParser

def build_area():
    area = ConfigParser.RawConfigParser()
    area.read('area.config')

    nx = int(area.get('Area', 'Nx'))
    ny = int(area.get('Area', 'Ny'))
    nz = int(area.get('Area', 'dNz'))
    R = ny/2 - 2
    x0 = nx/2
    y0 = ny/2
    z0 = nz/2

    r = float(area.get('Circle', 'R'))
    Y0 = 0.5
    Z0 = 0.5

    x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

    x_coord = np.loadtxt("area.x.coord", delimiter='\n')
    y_coord = np.loadtxt("area.y.coord", delimiter='\n')
    z_coord = np.loadtxt("area.z.coord", delimiter='\n')

    X = np.reshape(x_coord, (nx, 1, 1))
    Y = np.reshape(y_coord, (1, ny, 1))
    Z = np.reshape(z_coord, (1, 1, nz))

    in_boundary = (x == 0)
    out_boundary = (x == nx-1)

    mask = (x >= 0) & (x <= nx-1) &\
            (y >= 1) & (y <= ny-2) &\
            (z >= 1) & (z <= nz-3)

    inner_cylinder = ((Y-Y0)**2 + (Z-Z0)**2 < (r-0.1)**2)
    #inner_cylinder = ((y-y0)**2 + (z-z0)**2 < R**2)

    array = np.zeros((nx, ny, nz))

    array[inner_cylinder & mask] = 1
    array[in_boundary & inner_cylinder] = 2
    array[out_boundary & inner_cylinder] = 3

    """ Write mask file
    """
    with file('v_area.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
