#!/usr/bin/python
import sys
if len(sys.argv) < 2:
    sys.exit('Usage: %s file.vtk' % sys.argv[0])

from numpy import array
from time import sleep
from mayavi.api import Engine
from mayavi.modules.axes import Axes
from mayavi.modules.contour_grid_plane import ContourGridPlane
from mayavi.tools.show import show

try:
    engine = mayavi.engine
except NameError:
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()

scene = engine.scenes[0]
vtk_file_reader = engine.open(sys.argv[1], scene)
axes = Axes()
engine.add_filter(axes, vtk_file_reader)
contour_grid_plane = ContourGridPlane()
engine.add_filter(contour_grid_plane, vtk_file_reader)
contour_grid_plane.contour.filled_contours = True
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

contour_grid_plane_y = ContourGridPlane()
engine.add_filter(contour_grid_plane_y, vtk_file_reader)
contour_grid_plane_y.grid_plane.axis = 'y'
contour_grid_plane_y.grid_plane.position = 12
contour_grid_plane_y.actor.mapper.progress = 1.0
contour_grid_plane_y.actor.mapper.scalar_range = array([ 0.,  3.])
contour_grid_plane_y.actor.mapper.scalar_range = array([ 0.,  3.])
contour_grid_plane_y.actor.mapper.scalar_mode = 'use_cell_data'
contour_grid_plane_y.contour.filled_contours = True

contour_grid_plane.grid_plane.position = 30
show()
