#!/usr/bin/python
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
scene = engine.scenes[0]
vtk_file_reader = engine.open(u'/home/erthalion/coding/science/valves-model/build/surface.vtk', scene)
from mayavi.modules.axes import Axes
axes = Axes()
axes.property.color = (0.0, 0.0, 0.0)
engine.add_filter(axes, vtk_file_reader)
from mayavi.modules.surface import Surface
surface = Surface()
engine.add_filter(surface, vtk_file_reader)
scene.scene.background = (1.0, 1.0, 1.0)
scene.scene.foreground = (0.0, 0.0, 0.0)
surface.contour.contours = [1.5]
surface.actor.mapper.progress = 1.0
surface.actor.mapper.scalar_range = array([ 0.,  3.])
# surface.actor.mapper.input = <tvtk.tvtk_classes.poly_data.PolyData object at 0x8e826b0>
surface.enable_contours = True
# ------------------------------------------- 
from mayavi.tools.show import show
show()
