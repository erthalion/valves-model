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

    x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

    mask = (x >= 0) & (x <= nx-1) &\
            (y >= 0) & (y <= ny-2) &\
            (z >= 0) & (z <= ny-2)

    in_boundary = (x == 0)
    out_boundary = (x == nx-1)

    array = np.zeros((nx, ny, nz))

    array[mask] = 1
    array[in_boundary & mask] = 2
    array[out_boundary & mask] = 3

    """ Write mask file
    """
    with file('pressure.mask', 'w') as outfile:
        for slice_2d in array:
            np.savetxt(outfile, slice_2d, fmt="%i")
