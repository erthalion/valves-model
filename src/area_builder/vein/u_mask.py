#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import ConfigParser
from utils import take_off_first_last

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
    r_valves = 0.3
    X0 = 0.7
    Y0 = 0.5
    Z0 = 0.5

    x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

    x_coord = np.loadtxt("area.x.coord", delimiter='\n')
    y_coord = np.loadtxt("area.y.coord", delimiter='\n')
    z_coord = np.loadtxt("area.z.coord", delimiter='\n')

    X = np.reshape(x_coord, (nx, 1, 1))
    Y = np.reshape(y_coord, (1, ny, 1))
    Z = np.reshape(z_coord, (1, 1, nz))


    in_boundary = (x == 1)
    out_boundary = (x == nx-1)

    mask = (x >= 1) & (x <= nx-1) &\
            (y >= 1) & (y <= ny-3) &\
            (z >= 1) & (z <= ny-3)

    inner_cylinder = ((Y-Y0)**2 + (Z-Z0)**2 <= (r-0.01)**2)
    bottom_valve = ((X-X0)**2 + (Y-0*Y0)**2 <= (r_valves - 0.01)**2) &\
            (X < X0)

    top_valve = ((X-X0)**2 + (Y-2*Y0)**2 <= (r_valves - 0.01)**2) &\
            (X < X0)

    inner_cylinder = np.logical_and(
            np.roll(inner_cylinder, 1),
            np.roll(inner_cylinder, -1),
            )
    inner_cylinder = inner_cylinder.reshape((30, 30))
    inner_cylinder[inner_cylinder.any(axis=1).nonzero()[0][[0, -1]]] = False
    inner_cylinder = inner_cylinder.reshape((1, 30, 30))

#    bottom_valve = bottom_valve.reshape((30, 30))
    bottom_valve = np.logical_or(
            np.roll(bottom_valve, 1, axis=1),
            np.roll(bottom_valve, -1, axis=1),
            )

    bottom_valve = np.logical_or(
            np.roll(bottom_valve, 1, axis=0),
            bottom_valve
            )
    ##bottom_valve[bottom_valve.any(axis=1).nonzero()[0][[0, -1]]] = False
    ##bottom_valve = bottom_valve.reshape((1, 30, 30))

    ##top_valve = top_valve.reshape((30, 30))
    top_valve = np.logical_or(
            np.roll(top_valve, 1, axis=1),
            np.roll(top_valve, -1, axis=1),
            )
    top_valve = np.logical_or(
            np.roll(top_valve, 1, axis=0),
            top_valve
            )
    ##top_valve[top_valve.any(axis=1).nonzero()[0][[0, -1]]] = False
    ##top_valve = top_valve.reshape((1, 30, 30))

    array = np.zeros((nx, ny, nz))

    array[inner_cylinder & mask] = 1
    array[inner_cylinder & mask & bottom_valve] = 0
    array[inner_cylinder & mask & top_valve] = 0
    array[in_boundary & inner_cylinder] = 2
    array[out_boundary & inner_cylinder] = 3

    """ Write mask file
    """
    with file('u_area.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
