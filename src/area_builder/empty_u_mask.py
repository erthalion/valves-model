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

in_boundary = (x == 3)
out_boundary = (x == 58)

mask = (x > 1) & (x < 59) &\
        (y > 1) & (y < 24) &\
        (z > 1) & (z < 24)

array = np.zeros((nx, ny, nz))

array[mask] = 1
array[in_boundary & mask] = 2
array[out_boundary & mask] = 3

""" Write mask file
"""
with file('u_area.mask', 'w') as outfile:
    for slice_2d in array:
        np.savetxt(outfile, slice_2d, fmt="%i")
