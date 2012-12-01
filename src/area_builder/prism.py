import numpy as np
#import ipdb

#ipdb.set_trace()
nx = 25
ny = 25
nz = 61

x0 = 12
y0 = 12
r = 6
R = 12

x, y, z = np.ogrid[0:nx, 0:ny, 0:nz]

in_boundary = (0*x + 0*y + z == 1)
out_boundary = (0*x + 0*y + z == 60)
middle = (0*x +0*y + z > 1) & (0*x +0*y + z < 60)
middle_valve = (0*x +0*y + z > 20) & (0*x +0*y + z < 40)
valve = (0*x +0*y + z > 25) & (0*x +0*y + z < 30)
left = (0*x +0*y + z < 20)
right = (0*x +0*y + z > 40)
left_boundary = (0*x +0*y + z == 20)
right_boundary = (0*x + 0*y + z == 40)

inner_cylinder = ((x-x0)*(x-x0) + (y-y0)*(y-y0) + 0*z < r*r)

outer_prism = ((x-x0) + (y-y0) + 0*z <= R)\
        & (-(x-x0) + (y-y0) + 0*z <= R)\
        & ((x-x0) - (y-y0) + 0*z <= R)\
        & (-(x-x0) - (y-y0) + 0*z <= R)

mask = inner_cylinder & middle

left_mask = outer_prism & left

right_mask = outer_prism & right

first_valve = inner_cylinder\
        & (x + 0*y + 0*z < 12)\
        & (-(x-x0) + (y-y0) + 0*z > 0)\
        & middle_valve

second_valve = inner_cylinder\
        & (x + 0*y + 0*z > 12)\
        & ((x-x0) + (y-y0) + 0*z > 0)\
        & middle_valve

third_valve = inner_cylinder\
        & ((x-x0) + (y-y0) + 0*z < 0)\
        & (-(x-x0) + (y-y0) + 0*z < 0)\
        & middle_valve

input_mask = outer_prism & in_boundary
input_left_mask = outer_prism & left_boundary

output_mask = outer_prism & out_boundary
output_left_mask = outer_prism & right_boundary

array = np.zeros((nx, ny, nz))

array[mask] = 1
array[left_mask] = 1
array[right_mask] = 1
array[first_valve] = 0
array[second_valve] = 0
array[third_valve] = 0
array[input_mask] = 2
array[output_mask] = 3
array[input_left_mask] = 1
array[output_left_mask] = 1

output = open("mask.vtk","w")

header = "# vtk DataFile Version 1.0\n\
Data file for valves model\n\
ASCII\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS 61 25 25\n\
ORIGIN 0 0 0\n\
SPACING 1 1 1\n\
POINT_DATA 38125\n\
SCALARS scalars int\n\
LOOKUP_TABLE default\n"

output.write(header)

for i in range(0, 25):
    for j in range(0, 25):
        for k in range(0, 61):
            output.write(str(int(array.item((i, j, k))))+"\n")

output.close()
with file('prism.mask', 'w') as outfile:
    for slice_2d in array:
        np.savetxt(outfile, slice_2d, fmt="%i")
