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

    r = float(area.get('Circle', 'R'))
    r_valves = 0.3
    X0 = 0.7
    Y0 = 0.5
    Z0 = 0.5

    valve_width = 0.2
    width = 0.5
    height = 0.3

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
            (y >= 0) & (y <= ny-2) &\
            (z >= 0) & (z <= ny-2)

    inner_cylinder = ((Y-Y0)**2 + (Z-Z0)**2 <= (r-0.01)**2)

    bottom_valve = ((X-X0)**2 + (Y-0*Y0)**2 <= (r_valves - 0.01)**2) &\
            (X < X0)

    top_valve = ((X-X0)**2 + (Y-2*Y0)**2 <= (r_valves - 0.01)**2) &\
            (X < X0)

    array = np.zeros((nx, ny, nz))

    array[inner_cylinder & mask] = 1
    array[inner_cylinder & mask & bottom_valve] = 0
    array[inner_cylinder & mask & top_valve] = 0
    array[in_boundary & inner_cylinder] = 2
    array[out_boundary & inner_cylinder] = 3

    """ Write mask file
    """
    with file('pressure.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
