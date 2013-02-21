#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np

nx = 61
ny = 26
nz = 26

z0 = 13
y0 = 13
x0 = 30
r = 5
R = 11

x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

mask = (x > 0) & (x < 60) &\
        (y > 0) & (y < 25) &\
        (z > 0) & (z < 25)

array = np.zeros((nx, ny, nz))

array[mask] = 1

""" Write mask file
"""
with file('pressure.mask', 'w') as outfile:
    for slice_2d in array:
        np.savetxt(outfile, slice_2d, fmt="%i")