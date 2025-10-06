#%%
import os, sys
import codes
import numpy as np
import matplotlib.pyplot as plt

#%%
# MUST READ: 
# 1. fins out a s2 region that will initiate rotor by trails and errors
# 2. run main.py and set "rotor_flag = 0" to save the focal (s1 pacing) simulation's u and h values. 
#    NOTE that this step is necessary, because if "rotor_flag = 1", the u and h values will contain the s2 pacing stimulus, 
#    thus cannot find the correct u and h threshold for identifying a s2 pacing region.
# 3. run this python file to find out the u and h threshold.

# parameters for automatically find out s2 pacing regions. changes of these parameters will effect the rotor (location, shape, etc)
# the time interval after s1 pacing
s2_t = 60 

# the min max of the action potential and h determines the shape of the s2 pacing region
ap_min = (2.49 + 0.1) * 10**-7 # reference to np.min(action_potential_s2)
ap_max = (3.057 - 0.1) * 10**-6 # reference to np.max(action_potential_s2
h_min  = 0.015249441 + 0.003 # reference to np.min(h_s2)
h_max  = 0.027622167 - 0.003 # reference to np.max(h_s2)

#%%
# --------------------------------------------------
# load the rotor simulation figured out by manual trials and errors
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

folder_path = "../build/sim/bin/output/"
t, voltage, gate_variable, stress, xyz = codes.load_simulation_result.execute(folder_path)
# t[time_id]
# voltage[particles, time_steps]
# stress[particles, time_steps]
# xyz[particles, time_steps, coordinates]

#%%
# the manually assigned pacing sites
s1_pacing_particle_id = np.array([23403, 23409, 23410, 24111, 24112, 24113, 24118, 24119, 24120,
    24125, 24126, 24127, 24131, 24132, 24791, 24792, 24798, 24799, 24805, 24806, 24807])

n = 0 # time index
position = xyz[:,n,:]
s2_pacing_particle_id = np.array([], dtype=int)
for k in range(position.shape[0]):
    if position[k][0] >= 0.0 and position[k][0] <= 6.0:
        if position[k][1] >= -6.0:
            if position[k][2] >= 12.0:
                s2_pacing_particle_id = np.append(s2_pacing_particle_id,k)

# neighbor_id = neighbor_id_2d[s1_pacing_particle_id, :] # add all the neighbors of the pacing voxel to be paced
# neighbor_id = neighbor_id[neighbor_id != -1] # remove the -1s, which means no neighbors
# s1_pacing_particle_id = np.concatenate([s1_pacing_particle_id, neighbor_id])
# s1_pacing_particle_id = np.unique(s1_pacing_particle_id)

# neighbor_id = neighbor_id_2d[s2_pacing_particle_id, :] # add all the neighbors of the pacing voxel to be paced
# neighbor_id = neighbor_id[neighbor_id != -1] # remove the -1s, which means no neighbors
# s2_pacing_particle_id = np.concatenate([s2_pacing_particle_id, neighbor_id])
# s2_pacing_particle_id = np.unique(s2_pacing_particle_id)

# analyze the successful s2 region
action_potential_s2 = voltage[s2_pacing_particle_id,s2_t]
h_s2 = gate_variable[s2_pacing_particle_id,s2_t]

print(f"action potential min max: {np.min(action_potential_s2)} {np.max(action_potential_s2)}\n"
    f"h min max: {np.min(h_s2)} {np.max(h_s2)}")

# plot the values of action_potential_s2 and h_s2
plt.figure()
plt.plot(action_potential_s2,'b')
plt.plot(h_s2,'g')
plt.xlabel('particles')
plt.ylabel('')
plt.title('b: action potential, g: gate variable')

# plot the manually assigned pacing sites
codes.debug_display_of_s1s2_pacing_sites.execute(position, s1_pacing_particle_id, s2_pacing_particle_id)

#%%
# --------------------------------------------------
# automatically find s2 pacing voxels
action_potential_s2_t = voltage[:,s2_t]
h_s2_t = gate_variable[:,s2_t]

# automatically find s2 pacing sites
id1 = np.where((action_potential_s2_t >= ap_min) & (action_potential_s2_t <= ap_max))[0]
id2 = np.where((h_s2_t >= h_min) & (h_s2_t <= h_max))[0]
s2_pacing_particle_id_auto = np.intersect1d(id1, id2)

# plot the automatically assigned pacing sites
codes.debug_display_of_s1s2_pacing_sites.execute(position, s1_pacing_particle_id, s2_pacing_particle_id_auto)




# # grab a portion of the s2 pacing sites, so that it's like a curvy line instead of a ring
# id = s2_pacing_particle_id_auto[0] # find one voxel to start
# while id.size < s2_pacing_particle_id_auto.size/3: # repeat several times to include more neighbors
#     neighbor_id = neighbor_id_2d[id, :] # add all the neighbors of the pacing voxel to be paced
#     neighbor_id = neighbor_id[neighbor_id != -1] # remove the -1s, which means no neighbors
#     id = np.concatenate([np.atleast_1d(id), np.atleast_1d(neighbor_id)]) # add the neighbors
#     id = np.intersect1d(id, s2_pacing_particle_id_auto) # make sure its within the original shape

# # plot the automatically assigned pacing sites
# codes.debug_display_of_s1s2_pacing_sites.execute(voxel, s1_pacing_particle_id, id)

#%%
