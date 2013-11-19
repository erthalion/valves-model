import sys
import numpy as np
from pylab import quiver, plot, show

if len(sys.argv) != 2:
    print "File path is needed"
    sys.exit(1)

file_path = sys.argv[1]

x, y, x_vec, y_vec, x_ref, y_ref,\
    x1, y1, u1, v1,\
    x2, y2, u2, v2,\
    x3, y3, u3, v3,\
    x4, y4, u4, v4 = np.loadtxt(file_path, unpack=True)

grid_x = np.linspace(0, 2, 11)
grid_y = np.linspace(0, 1, 11)

# plot grid
X, Y = np.meshgrid(grid_x, grid_y)
plot(X, Y, 'k.', alpha=0.3)

# plot ref
plot(x_ref, y_ref, 'ro:')

# plot position
plot(x, y, 'ro-')

# plot vectors at points
quiver(x, y, x_vec, y_vec)

# plot point around
for i in range(x1.shape[0]):
    plot(
            x2[i], x3[i], x1[i], x2[i], x3[i], x4[i], x1[i],
            y2[i], y3[i], y1[i], y2[i], y3[i], y4[i], y1[i],
            ls=":"
            )

# plot vectors at the points

show()
