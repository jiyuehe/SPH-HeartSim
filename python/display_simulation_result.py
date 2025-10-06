# will need to install ffmpeg for saving the movie as a mp4 file
# sudo apt install ffmpeg

# %%
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import numpy as np

# %% load simulation results
folder_path = "../build/sim/bin/output/"
t, voltage, stress, xyz = fn_load_simulation_result.execute(folder_path)
# t[time_id]
# voltage[particles, time_steps]
# stress[particles, time_steps]
# xyz[particles, time_steps, coordinates]

debug_plot = 0
if debug_plot == 1:
    # plot the action potential voltage of a particle
    particle_id = 1800
    plt.figure()
    plt.plot(t, voltage[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('voltage')
    plt.title('Voltage Over Time of a Particle')
    plt.savefig('../result/voltage_of_a_particle.png')

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

# %% display simulation movie
movie_data = voltage

data_min = np.min(movie_data)
data_max = np.max(movie_data)
data_threshold = data_min
map_color = {}
n_time = movie_data.shape[1]
for n in range(n_time):
    if ((n+1) % (n_time//5)) == 0:
        print(f'compute color map {(n+1)/n_time*100:.1f}%')
    data = movie_data[:, n]
    color = convert_data_to_color.execute(data, data_min, data_max, data_threshold)
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
ax.view_init(elev = -50, azim = 100)
plot_handle = ax.scatter(voxel[:, 0], voxel[:, 1], voxel[:, 2], c=map_color[0], s=2, alpha=1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([x_min, x_max])
ax.set_ylim([y_min, y_max])
ax.set_zlim([z_min, z_max])
set_axes_equal.execute(ax)

pause_interval = 0.001
view_angles = {} # dictionary to store view angles for each frame
n_time = movie_data.shape[1]
for n in range(n_time):
    if ((n+1) % (n_time//5)) == 0:
        print(f'playing movie {(n+1)/n_time*100:.1f}%')

    plot_handle.set_facecolor(map_color[n]) # set color based on phase to each voxel
    ax.set_title(f'Time: {n}/{n_time} ms') # set title with current time step

    # capture current view angles
    elev = ax.elev  # elevation angle
    azim = ax.azim  # azimuth angle
    view_angles[n] = {'elev': elev, 'azim': azim}

    plt.pause(pause_interval)

# save simulation movie as mp4
if save_flag == 1:
    print("saving movie as mp4")

    def animate(n):
        if ((n+1) % (n_time//10)) == 0:
            print(f'saving movie {(n+1)/n_time*100:.1f}%')

        plot_handle.set_facecolor(map_color[n]) # set color based on phase to each voxel
        ax.set_title(f'Time: {n}/{n_time} ms') # set title with current time step

        ax.view_init(elev=view_angles[n]['elev'], azim=view_angles[n]['azim']) # restore view angle

    anim = animation.FuncAnimation(fig, animate, frames=n_time, interval=10, blit=False, repeat=False)
    # the interval parameter specifies the delay between frames in milliseconds

    # save
    writer = FFMpegWriter(fps=10, bitrate=1800)
    anim.save('result/simulation movie.mp4', writer=writer)

    print("movie saved as mp4")







n_time = xyz.shape[1]





do_flag = 1
if do_flag == 1:
    print("display movie")

    # dictionary to store view angles for each frame
    view_angles = {}

    interval = 0.01
    plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    for n in range(n_time):
        print(n/n_time)

        ax.clear()

        ax.scatter(xyz[:, n, 0], xyz[:, n, 1], xyz[:, n, 2], 
                   c=voltage[:, n], s=2, marker='.', alpha=1, cmap='coolwarm', vmin=v_min, vmax=v_max)
            
        # set title with current time step
        ax.set_title(f'Time: {t[n]}/{t[-1]}')
        
        # reset axis properties
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        fn_set_axes_equal.execute(ax)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])
        ax.set_zlim([z_min, z_max])

        # capture current view angles
        elev = ax.elev  # elevation angle
        azim = ax.azim  # azimuth angle
        view_angles[n] = {'elev': elev, 'azim': azim}

        plt.pause(interval)

    # save simulation movie as mp4
    do_flag = 1
    if do_flag == 1:
        print("saving movie as mp4")

        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection='3d')

        def animate(n):
            print(n/n_time)

            ax.clear()
            
            ax.scatter(xyz[:, n, 0], xyz[:, n, 1], xyz[:, n, 2], 
                    c=voltage[:, n], s=2, marker='.', alpha=1, cmap='coolwarm', vmin=v_min, vmax=v_max)
                
            # set title with current time step
            ax.set_title(f'Time: {t[n]}/{t[-1]}')
            
            # reset axis properties
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            fn_set_axes_equal.execute(ax)
            ax.set_xlim([x_min, x_max])
            ax.set_ylim([y_min, y_max])
            ax.set_zlim([z_min, z_max])

            # restore view angle to maintain user's rotation
            ax.view_init(elev=view_angles[n]['elev'], azim=view_angles[n]['azim'])

        anim = animation.FuncAnimation(fig, animate, frames=n_time, interval=10, blit=False, repeat=False)
        # the interval parameter specifies the delay between frames in milliseconds

        # save
        writer = FFMpegWriter(fps=10, bitrate=1800)
        anim.save('../result/simulation movie.mp4', writer=writer)

        print("movie saved as mp4")
