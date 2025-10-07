# will need to install ffmpeg for saving the movie as a mp4 file
# sudo apt install ffmpeg

# %%
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import numpy as np
import codes

# %% load simulation results
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

folder_path = "../build/sim/bin/output/"
t, voltage, gate_variable, stress, xyz = codes.load_simulation_result.execute(folder_path)
# t[time_id]
# voltage[particles, time_steps]
# stress[particles, time_steps]
# xyz[particles, time_steps, coordinates]

#%%
debug_plot = 0
if debug_plot == 1:
    particle_id = 1800

    # plot the action potential voltage of a particle
    plt.figure()
    plt.plot(t, voltage[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('voltage')
    plt.title('Voltage Over Time of a Particle')
    plt.savefig('../result/voltage_of_a_particle.png')

    # plot the gate variable of a particle
    plt.figure()
    plt.plot(t, gate_variable[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('gate variable')
    plt.title('Gate Variable Over Time of a Particle')
    plt.savefig('../result/gate_variable_of_a_particle.png')

    # plot the mechanical stress of a particle
    plt.figure()
    plt.plot(t, stress[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('stress')
    plt.title('Stress Over Time of a Particle')
    plt.savefig('../result/stress_of_a_particle.png')

    # plot a coordinate axis movements of a particle
    x_coordinate = [xyz[particle_id][time][0] for time in range(len(xyz[particle_id]))]
    y_coordinate = [xyz[particle_id][time][1] for time in range(len(xyz[particle_id]))]
    z_coordinate = [xyz[particle_id][time][2] for time in range(len(xyz[particle_id]))]
    fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    axs[0].plot(t, x_coordinate, 'r-', linewidth=1)
    axs[0].set_title('x coordinate movement')
    axs[1].plot(t, y_coordinate, 'g-', linewidth=1)
    axs[1].set_title('y coordinate movement')
    axs[2].plot(t, z_coordinate, 'b-', linewidth=1)
    axs[2].set_title('z coordinate movement')
    axs[2].set_xlabel('time, ms')
    plt.savefig('../result/movement_of_a_particle.png')

    plt.show()

# create phase from action potential
action_potential = voltage
v_gate = 0.13
action_potential_phase = np.zeros_like(action_potential)
activation_phase = np.zeros_like(action_potential)
for id in range(action_potential.shape[0]):
    if ((id+1) % (action_potential.shape[0]//5)) == 0:
        print(f'compute phase {(id+1)/action_potential.shape[0]*100:.1f}%')
    action_potential_phase[id,:], activation_phase[id,:] = codes.create_phase.execute(action_potential[id,:], v_gate)

debug_plot = 0
if debug_plot == 1:
    plt.figure()
    plt.plot(action_potential[1,:], 'b')
    plt.plot(action_potential_phase[1,:], 'g')
    plt.tight_layout()
    plt.show()

# %% display simulation movie
movie_data = action_potential_phase
data_min = np.min(movie_data)
data_max = np.max(movie_data)
data_threshold = data_min
map_color = {}
n_time = movie_data.shape[1]
for n in range(n_time):
    if ((n+1) % (n_time//5)) == 0:
        print(f'compute color map {(n+1)/n_time*100:.1f}%')
    data = movie_data[:, n]
    color = codes.convert_data_to_color.execute(data, data_min, data_max, data_threshold)
    map_color[n] = color

t_id = 0
d_buffer = 5
x_min = np.min(xyz[:,t_id,0]) - d_buffer
y_min = np.min(xyz[:,t_id,1]) - d_buffer
z_min = np.min(xyz[:,t_id,2]) - d_buffer
x_max = np.max(xyz[:,t_id,0]) + d_buffer
y_max = np.max(xyz[:,t_id,1]) + d_buffer
z_max = np.max(xyz[:,t_id,2]) + d_buffer

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection='3d')
ax.view_init(elev = 70, azim = -120)
n = 0 # time index
plot_handle = ax.scatter(xyz[:, n, 0], xyz[:, n, 1], xyz[:, n, 2], c=map_color[0], s=2, alpha=1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([x_min, x_max])
ax.set_ylim([y_min, y_max])
ax.set_zlim([z_min, z_max])
codes.set_axes_equal.execute(ax)

pause_interval = 0.001
view_angles = {} # dictionary to store view angles for each frame
n_time = movie_data.shape[1]
for n in range(n_time):
    if ((n+1) % (n_time//5)) == 0:
        print(f'playing movie {(n+1)/n_time*100:.1f}%')

    plot_handle._offsets3d = (xyz[:, n, 0], xyz[:, n, 1], xyz[:, n, 2]) # update voxels positions

    plot_handle.set_facecolor(map_color[n]) # set color based on phase to each voxel
    ax.set_title(f'Time: {round(t[n])}/{round(t[n_time-1])} ms') # set title with current time step

    # capture current view angles
    elev = ax.elev  # elevation angle
    azim = ax.azim  # azimuth angle
    view_angles[n] = {'elev': elev, 'azim': azim}

    plt.pause(pause_interval)

# save simulation movie as mp4
save_flag = 0
if save_flag == 1:
    print("saving movie as mp4")

    def animate(n):
        if ((n+1) % (n_time//10)) == 0:
            print(f'saving movie {(n+1)/n_time*100:.1f}%')

        plot_handle._offsets3d = (xyz[:, n, 0], xyz[:, n, 1], xyz[:, n, 2]) # update voxels positions

        plot_handle.set_facecolor(map_color[n]) # set color based on phase to each voxel
        ax.set_title(f'Time: {round(t[n])}/{round(t[n_time-1])} ms') # set title with current time step

        ax.view_init(elev=view_angles[n]['elev'], azim=view_angles[n]['azim']) # restore view angle

    anim = animation.FuncAnimation(fig, animate, frames=n_time, interval=10, blit=False, repeat=False)
    # the interval parameter specifies the delay between frames in milliseconds

    # save
    writer = FFMpegWriter(fps=10, bitrate=1800)
    anim.save('../result/simulation movie.mp4', writer=writer)

    print("movie saved as mp4")
