#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
#import ipdb

#ipdb.set_trace()
nx = 61
ny = 25
nz = 25

z0 = 12
y0 = 12
r = 5
R = 11

x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

in_fake_boundary = (0*z + 0*y + x == 0)
out_fake_boundary = (0*z + 0*y + x == nx -1)
in_boundary = (0*z + 0*y + x == 1)
out_boundary = (0*z + 0*y + x == nx - 2)
middle = (0*z +0*y + x > 1) & (0*z +0*y + x < nx - 1)
middle_valve = (0*z +0*y + x > 25) & (0*z +0*y + x < 28)
valve = (0*z +0*y + x > 25) & (0*z +0*y + x < 30)
left = (0*z +0*y + x < 20)
right = (0*z +0*y + x > 40)
left_boundary = (0*z +0*y + x == 20)
right_boundary = (0*z + 0*y + x == 40)

conuses_left = []
for i in range(6, 12):
    conus_mask = ((z-z0)*(z-z0) + (y-y0)*(y-y0) + 0*x < i*i) & (0*z +0*y + x > 25 + i) & (0*z +0*y + x < 48 + i)
    conuses_left.append(conus_mask)

conuses_right = []
for i in range(5, 12):
    conus_mask = ((z-z0)*(z-z0) + (y-y0)*(y-y0) + 0*x > i*i) & ((z-z0)*(z-z0) + (y-y0)*(y-y0) + 0*x < R*R) & (0*z +0*y + x > i - 5) & (0*z +0*y + x < 20 + i)
    conuses_right.append(conus_mask)

inner_cylinder = ((z-z0)*(z-z0) + (y-y0)*(y-y0) + 0*x < r*r)
outer_cylinder = ((z-z0)*(z-z0) + (y-y0)*(y-y0) + 0*x < R*R)

outer_prism = ((z-z0) + (y-y0) + 0*x < R)\
        & (-(z-z0) + (y-y0) + 0*x < R)\
        & ((z-z0) - (y-y0) + 0*x < R)\
        & (-(z-z0) - (y-y0) + 0*x < R)

mask = inner_cylinder & middle
outer_mask_left = outer_cylinder & left
outer_mask_right = outer_cylinder & right

left_mask = outer_prism & (0*z + 0*y + x < 20)

right_mask = outer_prism & right

first_valve = inner_cylinder\
        & (z + 0*y + 0*x < 12)\
        & (-(z-z0) + (y-y0) + 0*x > 0)\
        & middle_valve

second_valve = inner_cylinder\
        & (z + 0*y + 0*x > 12)\
        & ((z-z0) + (y-y0) + 0*x > 0)\
        & middle_valve

third_valve = inner_cylinder\
        & ((z-z0) + (y-y0) + 0*x < 0)\
        & (-(z-z0) + (y-y0) + 0*x < 0)\
        & middle_valve

input_mask = outer_cylinder & in_boundary
input_left_mask = outer_cylinder & left_boundary
input_fake_mask = outer_cylinder & in_fake_boundary

output_mask = outer_cylinder & out_boundary
output_left_mask = outer_cylinder & right_boundary
output_fake_mask = outer_cylinder & out_fake_boundary

array = np.zeros((nx, ny, nz))

array[mask] = 1
#array[outer_mask_left] = 1
#array[outer_mask_right] = 1
#array[left_mask] = 1
#array[right_mask] = 1
#array[first_valve] = 0
#array[second_valve] = 0
#array[third_valve] = 0

#array[input_left_mask] = 1
#array[output_left_mask] = 1
for conus in conuses_left:
    array[conus] = 1

for conus in conuses_right:
    array[conus] = 1

array[input_mask] = 2
array[output_mask] = 3
array[input_fake_mask] = 0

array[output_fake_mask] = 0

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
