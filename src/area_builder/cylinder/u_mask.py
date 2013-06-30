#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import ConfigParser
from math import tan, pi

def build_area():
    area = ConfigParser.RawConfigParser()
    area.read('area.config')

    nx = int(area.get('Area', 'Nx'))
    ny = int(area.get('Area', 'Ny'))
    nz = int(area.get('Area', 'dNz'))
    R = ny/2 - 2
    R_valves = ny/3
    x0 = nx/2
    y0 = ny/2
    z0 = nz/2
    alpha = pi/4
    big_valve_width = 6
    big_valve_shift = nx*2/3
    width = 1
    valves_open = 2

    x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

    in_boundary = (x == 1)
    out_boundary = (x == nx-1)

    mask = (x >= 1) & (x <= nx-1) &\
            (y >= 1) & (y <= ny-3) &\
            (z >= 1) & (z <= ny-3)

    inner_cylinder = ((y-y0)**2 + (z-z0)**2 < R**2)

    first_conuse = ((x-big_valve_shift) >= -np.sqrt((y-y0)**2 + (z-z0)**2)*tan(alpha)) &\
            (x-big_valve_shift <= -R/2) & (x-big_valve_shift >= -R)

    second_conuse = ((x-big_valve_shift - big_valve_width) >= -np.sqrt((y-y0)**2 + (z-z0)**2)*tan(alpha)) &\
            (x-big_valve_shift - big_valve_width <= -R/3) & (x-big_valve_shift - big_valve_width >= -R)

    first_valve = ((z-z0-valves_open) >= 0) & ((y-y0-valves_open) >= -0.5*(z-z0))
    second_valve = ((z-z0+valves_open) <= 0) & ((y-y0-valves_open) >= 0.5*(z-z0))
    third_valve = ((y-y0+valves_open) <= -0.5*(z-z0)) & ((y-y0+valves_open) <= 0.5*(z-z0))
    valves_width = ((x-nx*2/3+8)**2 + (y-y0)**2 + (z-z0)**2 >= (R_valves-width)**2) &\
            ((x-nx*2/3+8)**2 + (y-y0)**2 + (z-z0)**2 <= R_valves**2) &\
            (x-big_valve_shift >= -R/2)

    array = np.zeros((nx, ny, nz))

    array[inner_cylinder & mask] = 1
    array[in_boundary & inner_cylinder] = 2
    array[out_boundary & inner_cylinder] = 3

    """ Create big valve
    """
    array[inner_cylinder & first_conuse] = 0

    """ Create three valves
    """
    array[inner_cylinder & first_valve & valves_width] = 0
    array[inner_cylinder & second_valve & valves_width] = 0
    array[inner_cylinder & third_valve & valves_width] = 0

    """ Close big valve
    """
    array[inner_cylinder & second_conuse] = 1

    """ Write mask file
    """
    with file('u_area.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
