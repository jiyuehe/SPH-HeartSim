# %%
import matplotlib.pyplot as plt
import load_simulation_result

# %%
folder_path = "../build/sim/bin/output/"
t, voltage, stress, xyz = load_simulation_result.load(folder_path)

debug_plot = 0
if debug_plot == 1:
    # plot the action potential voltage of a particle
    particle_id = 1000
    plt.figure()
    plt.plot(t, voltage[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('voltage')
    plt.title('Voltage Over Time of a Particle')
    # plt.savefig('/home/j/Desktop/voltage_of_a_particle.png')

    # plot the mechanical stress of a particle
    particle_id = 1000
    plt.figure()
    plt.plot(t, stress[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('stress')
    plt.title('Stress Over Time of a Particle')
    # plt.savefig('/home/j/Desktop/stress_of_a_particle.png')

    # plot a coordinate axis movements of a particle
    particle_id = 1000
    x_coordinate = [xyz[particle_id][time][0] for time in range(len(xyz[particle_id]))]
    plt.figure()
    plt.plot(t, x_coordinate, 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('x coordinate')
    plt.title('x Coordinate vs Time of a Particle')
    # plt.savefig('/home/j/Desktop/x_coordinate_of_a_particle.png')


# %%
