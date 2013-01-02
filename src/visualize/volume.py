#!/usr/bin/python
import sys
if len(sys.argv) < 2:
    sys.exit('Usage: %s file.vtk' % sys.argv[0])

from numpy import array
from mayavi.modules.axes import Axes
from mayavi.api import Engine
from mayavi.modules.surface import Surface
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
axes.property.color = (0.0, 0.0, 0.0)
engine.add_filter(axes, vtk_file_reader)
surface = Surface()
engine.add_filter(surface, vtk_file_reader)
scene.scene.background = (1.0, 1.0, 1.0)
scene.scene.foreground = (0.0, 0.0, 0.0)
surface.contour.contours = [1.5]
surface.actor.mapper.progress = 1.0
surface.actor.mapper.scalar_range = array([ 0.,  3.])
surface.enable_contours = True

show()
