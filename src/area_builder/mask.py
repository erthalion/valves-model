#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np

nx = 61
ny = 25
nz = 25

z0 = 12
y0 = 12
x0 = 30
r = 5
R = 11

x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

in_boundary = (x == 2)
out_boundary = (x == 58)

mask = (x > 1) & (x < 59) &\
        (y > 1) & (y < 23) &\
        (z > 1) & (z < 23)

array = np.zeros((nx, ny, nz))

array[mask] = 1
array[in_boundary & mask] = 2
array[out_boundary & mask] = 3

""" Generate 3d array of coordinates
"""
Hx = Hy = Hz = 0.01
x_coord = np.zeros((nx))
y_coord = np.zeros((ny))
z_coord = np.zeros((nz))

for i in range(0, nx):
    x_coord[i] = i*Hx

for j in range(0, ny):
    y_coord[j] = j*Hy

for k in range(0, nz):
    z_coord[k] = k*Hz

output = open("mask.vtk","w")

header = "# vtk DataFile Version 1.0\n\
Data file for valves model\n\
ASCII\n\
DATASET STRUCTURED_GRID\n\
DIMENSIONS %d %d %d\n\
POINTS %d double\n" % (nx, ny, nz, nx*ny*nz)

output.write(header)

for k in range(0, nz):
    for j in range(0, ny):
        for i in range(0, nx):
            output.write('%0.10f %0.10f %0.10f\n' %(x_coord.item(i), y_coord.item(j), z_coord.item(k)))

data_header = "POINT_DATA %d\n\
SCALARS scalars int\n\
LOOKUP_TABLE default\n" % (nx*ny*nz)
output.write(data_header)
for k in range(0, nz):
    for j in range(0, ny):
        for i in range(0, nx):
            output.write('%d\n' % int(array.item((i, j, k))))

output.close()
