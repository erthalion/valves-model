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
    R = ny/2

    x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

    in_boundary = (x == 0)
    out_boundary = (x == nx-1)

    mask = (x >= 0) & (x <= nx-1) &\
            (y >= 1) & (y <= ny-3) &\
            (z >= 1) & (z <= nz-2)

    inner_cylinder = (y**2 + z**2 < R**2) & (x > 0)

    array = np.zeros((nx, ny, nz))

    array[in_boundary & mask & inner_cylinder] = 2
    array[out_boundary & mask & inner_cylinder] = 3
    array[inner_cylinder & mask] = 1

    """ Write mask file
    """
    with file('w_area.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
