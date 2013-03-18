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
    x0 = nx/2
    y0 = ny/2
    z0 = nz/2
    alpha = pi/4
    big_valve_width = 3
    big_valve_shift = nx*2/3

    x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

    in_boundary = (x == 0)
    out_boundary = (x == nx-1)

    mask = (x >= 0) & (x <= nx-1) &\
            (y >= 1) & (y <= ny-3) &\
            (z >= 1) & (z <= nz-2)

    inner_cylinder = ((y-y0)**2 + (z-z0)**2 < R**2)

    first_conuse = ((x-big_valve_shift) >= -np.sqrt((y-y0)**2 + (z-z0)**2)*tan(alpha)) &\
            (x-big_valve_shift <= -R/2) & (x-big_valve_shift >= -R)

    second_conuse = ((x-big_valve_shift - big_valve_width) >= -np.sqrt((y-y0)**2 + (z-z0)**2)*tan(alpha)) &\
            (x-big_valve_shift - big_valve_width <= -R/3) & (x-big_valve_shift - big_valve_width >= -R)

    array = np.zeros((nx, ny, nz))

    array[inner_cylinder & mask] = 1
    array[in_boundary & mask & inner_cylinder] = 2
    array[out_boundary & mask & inner_cylinder] = 3

    """ Create big valve
    """
    array[inner_cylinder & first_conuse] = 0
    array[inner_cylinder & second_conuse] = 1

    """ Write mask file
    """
    with file('w_area.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
