#!/usr/bin/python
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
# ------------------------------------------- 
scene = engine.scenes[0]
vtk_file_reader = engine.open(u'/home/erthalion/coding/science/valves-model/build/surface.vtk', scene)
axes = Axes()
engine.add_filter(axes, vtk_file_reader)
contour_grid_plane = ContourGridPlane()
engine.add_filter(contour_grid_plane, vtk_file_reader)
contour_grid_plane.contour.filled_contours = True
#scene.scene.camera.position = [0.0, 0.0, 0.03]
#scene.scene.camera.focal_point = [0.12, 0.0, 0.0]
#scene.scene.camera.view_angle = 30.0
#scene.scene.camera.view_up = [0.0, 0.3, 0.0]
#scene.scene.camera.clipping_range = [0.41681285720939437, 1.0957218138247224]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

contour_grid_plane_y = ContourGridPlane()
engine.add_filter(contour_grid_plane_y, vtk_file_reader)
contour_grid_plane_y.grid_plane.axis = 'y'
contour_grid_plane_y.grid_plane.position = 12
contour_grid_plane_y.actor.mapper.progress = 1.0
contour_grid_plane_y.actor.mapper.scalar_range = array([ 0.,  3.])
# contour_grid_plane.actor.mapper.input = <tvtk.tvtk_classes.poly_data.PolyData object at 0x19497f50>
contour_grid_plane_y.actor.mapper.scalar_range = array([ 0.,  3.])
contour_grid_plane_y.actor.mapper.scalar_mode = 'use_cell_data'
contour_grid_plane_y.contour.filled_contours = True

contour_grid_plane.grid_plane.position = 30
#for i in range(0, 60):
    #contour_grid_plane.grid_plane.position = i
    #sleep(0.2)
# ------------------------------------------- 
show()
