# %%
import os
import re
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import random

script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

# %% find all the simulation output files
folder_path="../build/sim/bin/output/"
all_files = os.listdir(folder_path)
file_names = [f for f in all_files if f.startswith("PhysiologyHeart_")]

# sort the filenames
numbers = []
i = 0
for file_name in file_names:
    match = re.search(r'PhysiologyHeart_(\d+)', file_name)
    numbers.append(int(match.group(1)))

paired = list(zip(numbers, file_names)) # create pairs 
paired.sort()
numbers_sorted = [pair[0] for pair in paired]
file_names_sorted = [pair[1] for pair in paired]

# load the electro-mechanical simulation result
xyz_temp = []
voltage_temp = []
stress_temp = []
for file_name in file_names_sorted:
    tree = ET.parse(folder_path + file_name)
    root = tree.getroot()

    # load xyz position of the particles
    data_array = None
    for da in root.iter("DataArray"):
        if da.attrib.get("Name") == "Position":
            data_array = da
            break

    text_data = data_array.text.strip().split() # extract the ASCII text values
    data = np.array(list(map(float, text_data)), dtype=float)
    xyz_temp.append(data.reshape(-1, 3)) # reshape to N x 3 (x,y,z)

    # load voltage of the particles
    data_array = None
    for da in root.iter("DataArray"):
        if da.attrib.get("Name") == "Voltage":
            data_array = da
            break

    text_data = data_array.text.strip().split() # extract the ASCII text values
    data = np.array(list(map(float, text_data)), dtype=float)
    voltage_temp.append(data)

    # load stress of the particles
    data_array = None
    for da in root.iter("DataArray"):
        if da.attrib.get("Name") == "ActiveContractionStress":
            data_array = da
            break

    text_data = data_array.text.strip().split() # extract the ASCII text values
    data = np.array(list(map(float, text_data)), dtype=float)
    stress_temp.append(data)

# save data in time series
t = [] # time 
for n in numbers_sorted: # n: time stamps
    t.append(n/1000000) # time

voltage = []
num_particles = len(voltage_temp[0])
num_time_steps = len(voltage_temp)
for particle in range(num_particles):
    particle_voltages = []
    for time in range(num_time_steps):
        particle_voltages.append(voltage_temp[time][particle])
    voltage.append(particle_voltages)

debug_plot = 0
if debug_plot == 1:
    particle_id = 1000
    plt.figure()
    plt.plot(t, voltage[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('voltage')
    plt.title('Voltage Over Time for a Particle')
    plt.savefig('/home/j/Desktop/voltage_of_a_particle.png')

stress = []
num_particles = len(stress_temp[0])
num_time_steps = len(stress_temp)
for particle in range(num_particles):
    particle_stress = []
    for time in range(num_time_steps):
        particle_stress.append(stress_temp[time][particle])
    stress.append(particle_stress)

debug_plot = 0
if debug_plot == 1:
    particle_id = 1000
    plt.figure()
    plt.plot(t, stress[particle_id], 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('stress')
    plt.title('Stress Over Time for a Particle')
    plt.savefig('/home/j/Desktop/stress_of_a_particle.png')

xyz = [] # xyz[particle_id][time][coordinate_id]
num_particles = len(xyz_temp[0])
num_time_steps = len(xyz_temp)
for particle in range(num_particles):
    particle_xyz = []
    for time in range(num_time_steps):
        particle_xyz.append(xyz_temp[time][particle])
    xyz.append(particle_xyz)

debug_plot = 0
if debug_plot == 1:
    particle_id = 1000
    x_coordinate = [xyz[particle_id][time][0] for time in range(len(xyz[particle_id]))]
    plt.figure()
    plt.plot(t, x_coordinate, 'b-', linewidth=1)
    plt.xlabel('time, unit: ms')
    plt.ylabel('x coordinate')
    plt.title('x Coordinate vs Time for a Particle')
    plt.savefig('/home/j/Desktop/x_coordinate_of_a_particle.png')
