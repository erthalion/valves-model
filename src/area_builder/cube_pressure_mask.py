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

in_boundary = (0*z + 0*y + x == 11)
out_boundary = (0*z + 0*y + x == 49)

inner_cube = (0*(z-z0) + 0*(y-y0) + np.fabs(x-x0) < 5) &\
                (0*(z-z0) + np.fabs(y-y0) + 0*(x-x0) < 5) &\
                (np.fabs(z-z0) + 0*(y-y0) + 0*(x-x0) < 5)

outer_cube = (0*(z-z0) + 0*(y-y0) + np.fabs(x-x0) < 20) &\
                (0*(z-z0) + np.fabs(y-y0) + 0*(x-x0) < 10) &\
                (np.fabs(z-z0) + 0*(y-y0) + 0*(x-x0) < 10)

input_mask = outer_cube & in_boundary

output_mask = outer_cube & out_boundary

array = np.zeros((nx, ny, nz))

array[outer_cube] = 1
array[inner_cube] = 0
array[input_mask] = 2
array[output_mask] = 3

""" Write mask file
"""
with file('pressure.mask', 'w') as outfile:
    for slice_2d in array:
        np.savetxt(outfile, slice_2d, fmt="%i")

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

output = open("pressure_mask.vtk","w")

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
