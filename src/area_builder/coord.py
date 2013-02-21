#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np

nx = 61
ny = 26
nz = 26

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
