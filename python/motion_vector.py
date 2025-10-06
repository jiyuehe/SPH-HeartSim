# %%
import os
import matplotlib.pyplot as plt
import numpy as np
import codes

# %%
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

folder_path = "../build/sim/bin/output/"
t, voltage, stress, xyz = codes.load_simulation_result.execute(folder_path)
# t[time_id]
# voltage[particles, time_steps]
# stress[particles, time_steps]
# xyz[particles, time_steps, coordinates]

particles_1 = xyz[:, 50, :]
particles_2 = xyz[:, 52, :]

debug_plot = 0
if debug_plot == 1:
    # display together the particles location of two different times
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(particles_1[:, 0], particles_1[:, 1], particles_1[:, 2], color='blue', s=1, marker='.')
    ax.scatter(particles_2[:, 0], particles_2[:, 1], particles_2[:, 2], color='red', s=1, marker='.')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,1])
    codes.set_axes_equal.execute(ax)
    plt.show()

# %%
# compute the motion vector
motion_vector = particles_2 - particles_1

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
    codes.set_axes_equal.execute(ax)
    plt.show()

# %%
