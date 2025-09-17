# to activate a Python virtual environment
# source /home/j/myenv/bin/activate
# to deactivate the virtual environment
# deactivate

# %%
import os
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import fn_set_axes_equal

# %% load the xml file
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

tree = ET.parse('../build/sim/bin/reload/HeartModel_rld.xml')
root = tree.getroot()

N = len(root)
OriginalID = np.zeros(N) # 1xN
HeartModelPart2ID = np.zeros(N)
VolumetricMeasure = np.zeros(N)
Phi = np.zeros(N)
Position = np.zeros((N, 3)) # Nx3
Fiber = np.zeros((N, 3))
Sheet = np.zeros((N, 3))
i = 0
for child in root:
    # print(child.attrib)
    OriginalID[i] = child.attrib['OriginalID']
    HeartModelPart2ID[i] = child.attrib['HeartModelPart2ID']
    VolumetricMeasure[i] = child.attrib['VolumetricMeasure']
    Phi[i] = child.attrib['Phi']
    
    temp = child.attrib['Position']
    Position[i,:] = np.array([float(x.strip()) for x in temp.split(',')])

    temp = child.attrib['Fiber']
    Fiber[i,:] = np.array([float(x.strip()) for x in temp.split(',')])

    temp = child.attrib['Sheet']
    Sheet[i,:] = np.array([float(x.strip()) for x in temp.split(',')])

    i = i+1

x_min, y_min, z_min = np.min(Position, axis=0)
x_max, y_max, z_max = np.max(Position, axis=0)
print(f"X range: {x_min:.3f} to {x_max:.3f}")
print(f"Y range: {y_min:.3f} to {y_max:.3f}")
print(f"Z range: {z_min:.3f} to {z_max:.3f}")

debug_plot = 0
if debug_plot == 1:
    # plot the particles of HeartModelPart2ID == 0 and HeartModelPart2ID == 2
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    id1 = HeartModelPart2ID == 0
    id2 = HeartModelPart2ID == 2
    ax.scatter(Position[id1, 0], Position[id1, 1], Position[id1, 2], color='blue', label='ID==0', s=1, marker='.')
    ax.scatter(Position[id2, 0], Position[id2, 1], Position[id2, 2], color='red', label='ID==2', s=1, marker='.')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1,1,1])
    fn_set_axes_equal.execute(ax)
    plt.legend()

# %% pacing site
pacing_id_s1 = np.zeros(0)
i = 0
for p in Position:
    # print(p)
    if -4.0 <= p[0] and p[0] <= 0.0: # x coordinates
        if -20.0 <= p[1] and p[1] <= -12.0: # y coordinates
            if 136.0 <= p[2] and p[2] <= 140.0: # z coordinates
                pacing_id_s1 = np.append(pacing_id_s1,i)
    i = i+1
pacing_id_s1 = pacing_id_s1.astype(int)

pacing_id_s2 = np.zeros(0)
i = 0
for p in Position:
    # print(p)
    if 0.0 <= p[0] and p[0] <= 6.0: # x coordinates
        if -6.0 <= p[1]: #and p[1] <= -38.0: # y coordinates
            if 12.0 <= p[2]: #and p[2] <= 2.0: # z coordinates
                pacing_id_s2 = np.append(pacing_id_s2,i)
    i = i+1

pacing_id_s2 = pacing_id_s2.astype(int)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Position[:, 0], Position[:, 1], Position[:, 2], color='blue', s=1, marker='.')
ax.scatter(Position[pacing_id_s1, 0], Position[pacing_id_s1, 1], Position[pacing_id_s1, 2], color='red', s=50, marker='.')
ax.scatter(Position[pacing_id_s2, 0], Position[pacing_id_s2, 1], Position[pacing_id_s2, 2], color='g', s=50, marker='.')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_box_aspect([1,1,1])
fn_set_axes_equal.execute(ax)
plt.show()
