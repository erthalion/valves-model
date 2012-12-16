# Recorded script from Mayavi2
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
vtk_file_reader = engine.open(u'/home/erthalion/coding/cpp/vtk/data.vtk')
from mayavi.modules.axes import Axes
axes = Axes()
engine.add_filter(axes, vtk_file_reader)
from mayavi.modules.streamline import Streamline
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
#streamline1.seed.widget.handle_direction = array([ 1.,  0.,  0.])
streamline1.seed.widget.enabled = False
streamline1.seed.widget.center = array([ 30.,  12.,  12.])
#streamline1.seed.widget.handle_direction = array([ 1.,  0.,  0.])
streamline1.seed.widget.interactor = None
# streamline1.seed.widget = <tvtk.tvtk_classes.plane_widget.PlaneWidget object at 0xd0f2410>
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
from mayavi.tools.show import show
show()
