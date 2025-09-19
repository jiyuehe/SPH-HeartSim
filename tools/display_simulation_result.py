# will need to install ffmpeg for saving the movie as a mp4 file
# sudo apt install ffmpeg

# %%
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import numpy as np
import fn_load_simulation_result # function of loading simulation results
import fn_set_axes_equal

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
    particle_id = 500
    plt.figure()
    plt.plot(t, voltage[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('voltage')
    plt.title('Voltage Over Time of a Particle')
    # plt.savefig('/home/j/Desktop/voltage_of_a_particle.png')

    # plot the mechanical stress of a particle
    plt.figure()
    plt.plot(t, stress[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('stress')
    plt.title('Stress Over Time of a Particle')
    # plt.savefig('/home/j/Desktop/stress_of_a_particle.png')

    # plot a coordinate axis movements of a particle
    x_coordinate = [xyz[particle_id][time][0] for time in range(len(xyz[particle_id]))]
    plt.figure()
    plt.plot(t, x_coordinate, 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('x coordinate')
    plt.title('x Coordinate vs Time of a Particle')
    # plt.savefig('/home/j/Desktop/x_coordinate_of_a_particle.png')

    plt.show()

# %% display simulation movie
num_particles, num_time_steps, _ = xyz.shape
v_min = np.min(voltage)
v_max = np.max(voltage)

t_id = 0
d_buffer = 5
x_min = np.min(xyz[:,t_id,0]) - d_buffer
y_min = np.min(xyz[:,t_id,1]) - d_buffer
z_min = np.min(xyz[:,t_id,2]) - d_buffer
x_max = np.max(xyz[:,t_id,0]) + d_buffer
y_max = np.max(xyz[:,t_id,1]) + d_buffer
z_max = np.max(xyz[:,t_id,2]) + d_buffer

do_flag = 1
if do_flag == 1:
    print("display movie")

    # dictionary to store view angles for each frame
    view_angles = {}

    interval = 0.01
    plt.figure(figsize=(10, 8))
    ax = plt.axes(projection='3d')
    for n in range(num_time_steps):
        print(n/num_time_steps)

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
            print(n/num_time_steps)

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

        anim = animation.FuncAnimation(fig, animate, frames=num_time_steps, interval=10, blit=False, repeat=False)
        # the interval parameter specifies the delay between frames in milliseconds

        # save
        writer = FFMpegWriter(fps=10, bitrate=1800)
        anim.save('simulation movie.mp4', writer=writer)

        print("movie saved as mp4")
