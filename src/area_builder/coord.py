#!/usr/bin/env python
# -*- coding: utf8 -*-

import numpy as np
import ConfigParser

def build_area():
    area = ConfigParser.RawConfigParser()
    area.read('area.config')

    nx = int(area.get('Area', 'Nx'))
    ny = int(area.get('Area', 'Ny'))
    nz = int(area.get('Area', 'dNz'))

    lx = float(area.get('Area', 'Lx'))
    ly = float(area.get('Area', 'Ly'))
    lz = float(area.get('Area', 'Lz'))

    """ Generate 3d array of coordinates
    """
    Hx = lx / nx
    Hy = ly / ny
    Hz = lz / nz
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
