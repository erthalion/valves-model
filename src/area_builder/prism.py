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
out_boundary = (0*x + 0*y + z == 61)
middle = (0*x +0*y + z > 1) & (0*x +0*y + z < 60)
left = (0*x +0*y + z < 20)
right = (0*x +0*y + z > 40)
not_left_boundary = (0*x +0*y + z != 20)
not_right_boundary = (0*x + 0*y + z != 40)

prism = ((x-x0)*(x-x0) + (y-y0)*(y-y0) + 0*z < r*r)
second_prism = ((x-x0) + (y-y0) + 0*z < R) & (-(x-x0) + (y-y0) + 0*z < R) & ((x-x0) - (y-y0) + 0*z < R) & (-(x-x0) - (y-y0) + 0*z < R)
mask = prism & middle
left_mask = second_prism & left
right_mask = second_prism & right

first_valve = prism & (x + 0*y + 0*z < 12) & (-(x-x0) + (y-y0) + 0*z > 0) & not_left_boundary & not_right_boundary & middle
second_valve = prism & (x + 0*y + 0*z > 12) & ((x-x0) + (y-y0) + 0*z > 0) & not_left_boundary & not_right_boundary & middle
third_valve = prism & ((x-x0) + (y-y0) + 0*z < 0) & (-(x-x0) + (y-y0) + 0*z < 0) & not_left_boundary & not_right_boundary & middle

input_mask = second_prism & in_boundary
output_mask = second_prism & out_boundary

array = np.zeros((nx, ny, nz))

array[mask] = 1
array[left_mask] = 1
array[right_mask] = 1
array[first_valve] = 0
array[second_valve] = 0
array[third_valve] = 0
array[input_mask] = 2
array[output_mask] = 2

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
