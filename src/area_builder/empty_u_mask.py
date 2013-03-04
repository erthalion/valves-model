#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import ConfigParser

area = ConfigParser.RawConfigParser()
area.read('area.config')

nx = int(area.get('Area', 'Nx'))
ny = int(area.get('Area', 'Ny'))
nz = int(area.get('Area', 'dNz'))

x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

in_boundary = (x == 1)
out_boundary = (x == nx-1)

mask = (x >= 1) & (x <= nx-1) &\
        (y >= 1) & (y <= ny-3) &\
        (z >= 1) & (z <= ny-3)

array = np.zeros((nx, ny, nz))

array[mask] = 1
array[in_boundary & mask] = 2
array[out_boundary & mask] = 3

""" Write mask file
"""
with file('u_area.mask', 'w') as outfile:
    for slice_2d in array:
        np.savetxt(outfile, slice_2d, fmt="%i")
