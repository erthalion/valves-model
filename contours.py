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
vtk_file_reader = engine.open(u'/home/erthalion/coding/science/valves-model/build/surface.vtk', scene)
from mayavi.modules.axes import Axes
axes = Axes()
engine.add_filter(axes, vtk_file_reader)
from mayavi.modules.contour_grid_plane import ContourGridPlane
contour_grid_plane = ContourGridPlane()
engine.add_filter(contour_grid_plane, vtk_file_reader)
scene.scene.camera.position = [-0.57811328100914205, 0.14564046162201691, -0.070405212439518763]
scene.scene.camera.focal_point = [0.12, 0.0, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.090547904527663176, -0.046866218917022473, -0.99478873863251538]
scene.scene.camera.clipping_range = [0.41681285720939437, 1.0957218138247224]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
contour_grid_plane.grid_plane.position = 16
contour_grid_plane.grid_plane.position = 10
contour_grid_plane.grid_plane.position = 4
contour_grid_plane.grid_plane.position = 0
contour_grid_plane.grid_plane.position = 6
contour_grid_plane.grid_plane.position = 12
contour_grid_plane.grid_plane.position = 18
contour_grid_plane.grid_plane.position = 24
contour_grid_plane.grid_plane.position = 30
contour_grid_plane.grid_plane.position = 36
contour_grid_plane.grid_plane.position = 42
contour_grid_plane.grid_plane.position = 48
contour_grid_plane.grid_plane.position = 54
contour_grid_plane.grid_plane.position = 60
# ------------------------------------------- 
from mayavi.tools.show import show
show()
