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

in_boundary = (0*z + 0*y + x == 1)
out_boundary = (0*z + 0*y + x == 59)

mask = (0*z + 0*y + x > 0) & (0*z + 0*y + x < 60) &\
        (0*z + y + 0*x > 0) & (0*z + y + 0*x < 24) &\
        (z + 0*y + 0*x > 0) & (z + 0*y + 0*x < 24)

array = np.zeros((nx, ny, nz))

array[mask] = 1
array[in_boundary & mask] = 2
array[out_boundary & mask] = 3

""" Write mask file
"""
with file('area.mask', 'w') as outfile:
    for slice_2d in array:
        np.savetxt(outfile, slice_2d, fmt="%i")

""" Generate 3d array of coordinates
"""
Hx = Hy = Hz = 0.01
x_coord = np.zeros((nx))
y_coord = np.zeros((ny))
z_coord = np.zeros((nz))

output = open("area.x.coord","w")
for i in range(0, nx):
    x_coord[i] = i*Hx
    output.write('%0.10f\n' % x_coord.item((i)))
output.close()

output = open("area.y.coord","w")
for j in range(0, ny):
    y_coord[j] = j*Hy
    output.write('%0.10f\n' % y_coord.item((j)))
output.close()

output = open("area.z.coord","w")
for k in range(0, nz):
    z_coord[k] = k*Hz
    output.write('%0.10f\n' % z_coord.item((k)))
output.close()

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
