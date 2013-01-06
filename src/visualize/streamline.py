#!/usr/bin/python
import sys
if len(sys.argv) < 2:
    sys.exit('Usage: %s file.vtk' % sys.argv[0])

from numpy import array
from mayavi.modules.axes import Axes
from mayavi.api import Engine
from mayavi.modules.streamline import Streamline
from mayavi.tools.show import show

try:
    engine = mayavi.engine
except NameError:
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()

scene = engine.scenes[0]
vtk_file_reader = engine.open(sys.argv[1])
axes = Axes()
engine.add_filter(axes, vtk_file_reader)
streamline = Streamline()
engine.add_filter(streamline, vtk_file_reader)
streamline.seed.widget.center = array([  0.,  12.,  12.])
streamline.seed.widget.center = array([  0.,  12.,  12.])
streamline.seed.widget.handle_direction = array([ 1.,  0.,  0.])
streamline.seed.widget.enabled = False

streamline1 = Streamline()
engine.add_filter(streamline1, vtk_file_reader)
streamline1.seed.widget = streamline1.seed.widget_list[2]
streamline1.seed.widget.center = array([ 30.,  12.,  12.])
streamline1.seed.widget.enabled = False
streamline1.seed.widget.center = array([ 30.,  12.,  12.])
streamline1.seed.widget.interactor = None
streamline1.seed.widget.center = array([ 39.,  12.,  12.])
streamline1.seed.widget.origin = array([ 39.,   6.,   6.])
streamline1.seed.widget.origin = array([ 39.,   6.,   6.])
streamline1.seed.widget.center = array([ 39.,  12.,  12.])
streamline1.seed.widget.normal = array([ 1.,  0.,  0.])
streamline1.seed.widget.point1 = array([ 39.,  18.,   6.])
streamline1.seed.widget.point1 = array([ 39.,  18.,   6.])
streamline1.seed.widget.point2 = array([ 39.,   6.,  18.])
streamline1.seed.widget.point2 = array([ 39.,   6.,  18.])
streamline1.seed.widget.enabled = False
scene.scene.background = (1.0, 1.0, 1.0)
scene.scene.foreground = (0.0, 0.0, 0.0)
show()
