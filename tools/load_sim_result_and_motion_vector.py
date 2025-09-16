# %%
import os
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import random

script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

# %% load the electro-mechanical simulation result
tree = ET.parse('../build/sim/bin/output/PhysiologyHeart_0028654166.vtp')
root = tree.getroot()

# find the DataArray with Name="Position"
data_array = None
for da in root.iter("DataArray"):
    if da.attrib.get("Name") == "Position":
        data_array = da
        break

text_data = data_array.text.strip().split() # extract the ASCII text values
coords = np.array(list(map(float, text_data)), dtype=float)
particles_1 = coords.reshape(-1, 3) # reshape to N x 3 (x,y,z)

tree = ET.parse('../build/sim/bin/output/PhysiologyHeart_0026449999.vtp')
root = tree.getroot()

# find the DataArray with Name="Position"
data_array = None
for da in root.iter("DataArray"):
    if da.attrib.get("Name") == "Position":
        data_array = da
        break

text_data = data_array.text.strip().split() # extract the ASCII text values
coords = np.array(list(map(float, text_data)), dtype=float)
particles_2 = coords.reshape(-1, 3) # reshape to N x 3 (x,y,z)

# plot
debug_plot = 0
if debug_plot == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(particles_1[:, 0], particles_1[:, 1], particles_1[:, 2], color='blue', s=1, marker='.')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,1])
    plt.legend()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(particles_2[:, 0], particles_2[:, 1], particles_2[:, 2], color='blue', s=1, marker='.')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,1])
    plt.legend()

# %%
# compute the motion vector
motion_vector = particles_2 - particles_1

def set_axes_equal(ax):
    """Make 3D plot axes have equal scale so spheres look like spheres."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max([x_range, y_range, z_range]) / 2.0

    mid_x = np.mean(x_limits)
    mid_y = np.mean(y_limits)
    mid_z = np.mean(z_limits)

    ax.set_xlim3d([mid_x - max_range, mid_x + max_range])
    ax.set_ylim3d([mid_y - max_range, mid_y + max_range])
    ax.set_zlim3d([mid_z - max_range, mid_z + max_range])

# plot
debug_plot = 1
if debug_plot == 1:
    N = 1500
    step = len(particles_1) // N # at most ~N arrows
    X, Y, Z = particles_1[::step].T # array[start:stop:step]
    U, V, W = motion_vector[::step].T
    # ids = random.sample(range(0, len(particles_1)), N)
    # X, Y, Z = particles_1[ids].T 
    # U, V, W = motion_vector[ids].T

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.quiver(X, Y, Z, U, V, W, length = 40.0, normalize = False, color="blue")
    # ax.scatter(particles_1[:, 0], particles_1[:, 1], particles_1[:, 2], color='gray', s=15, marker='.')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Motion vectors")
    ax.set_box_aspect([1,1,1])  # equal scale for x, y, z (this does not work)
    set_axes_equal(ax)   # set equal aspect ratio
    plt.show()
