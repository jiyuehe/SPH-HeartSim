# %%
import os
import re
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

def load(folder_path):
    script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
    os.chdir(script_dir) # change the working directory

    # %% find all the simulation output files
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

    # %% save data in time series
    t = np.array(numbers_sorted) / 1000000  # time in ms
    voltage = np.array(voltage_temp).T # .T is transpose
    stress = np.array(stress_temp).T
    xyz = np.array(xyz_temp).transpose(1, 0, 2) # xyz[particles, time_steps, coordinates]
        # .transpose(1, 0, 2) means:
        # xyz_temp's axis 0 (time axis) is moved to position 1
        # xyz_temp's axis 1 (particle axis) is moved to position 0
        # xyz_temp's axis 2 (coordinate axis) stays in position 2

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

    return t, voltage, stress, xyz
