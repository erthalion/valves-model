#!/usr/bin/python
import sys
if len(sys.argv) < 3:
    sys.exit('Usage: %s mask.vtk streamlines.vtk' % sys.argv[0])

from numpy import array
from mayavi.api import Engine
from mayavi.modules.axes import Axes
from mayavi.modules.surface import Surface
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
vtk_file_reader = engine.open(sys.argv[1], scene)
vtk_file_reader1 = engine.open(sys.argv[2], scene)

surface = Surface()
engine.add_filter(surface, vtk_file_reader)
surface.contour.contours = [1.5]
surface.actor.mapper.progress = 1.0
surface.actor.mapper.scalar_range = array([ 0.,  3.])
surface.enable_contours = True
surface.actor.property.point_size = 0.0
surface.actor.property.opacity = 0.1602

streamline = Streamline()
streamline.seed.widget = streamline.seed.widget_list[2]
engine.add_filter(streamline, vtk_file_reader1)
streamline.stream_tracer.integration_direction = 'both'
#streamline.seed.widget.enabled = False
#streamline.seed.widget.interactor = None
streamline.actor.property.specular_color = (0.0, 0.0, 0.0)
streamline.actor.property.diffuse_color = (0.0, 0.0, 0.0)
streamline.actor.property.ambient_color = (0.0, 0.0, 0.0)
streamline.actor.property.color = (0.0, 0.0, 0.0)
#streamline.seed.widget.origin = array([ 0.00733424, -0.01543501, -0.03262647])
#streamline.seed.widget.center = array([ 0.00733424,  0.12744006,  0.11075326])
#streamline.seed.widget.normal = array([ 1.,  0.,  0.])
#streamline.seed.widget.point1 = array([ 0.00733424,  0.27031514, -0.03262647])
#streamline.seed.widget.point2 = array([ 0.00733424, -0.01543501,  0.25413299])
#streamline.seed.widget.enabled = False
streamline.seed.widget.resolution = 4

axes = Axes()
engine.add_filter(axes, vtk_file_reader)
#import ipdb; ipdb.set_trace()
axes.property.reference_count = 17
axes.axes.property.color = (0.0, 0.0, 0.0)
axes.title_text_property.shadow_offset = array([ 1, -1])
axes.title_text_property.color = (0.0, 0.0, 0.0)
axes.label_text_property.shadow_offset = array([ 1, -1])
axes.label_text_property.color = (0.0, 0.0, 0.0)


scene.scene.background = (1.0, 1.0, 1.0)
scene.scene.disable_render = True
scene.scene.disable_render = False
scene.scene.camera.position = [0.32692494007520589, 1.5773625791214541, 0.013388461834330156]
scene.scene.camera.focal_point = [0.29999999999999999, 0.12199067883193493, 0.11019065231084824]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.999824689874523, -0.018625038896002939, -0.0019228737443115677]
scene.scene.camera.clipping_range = [1.0908954355408496, 1.92551655925284]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
#scene.scene.save(u'/home/erthalion/snapshot.png')

show()
