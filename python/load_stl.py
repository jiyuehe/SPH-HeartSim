# %%
import os
from stl import mesh # pip install numpy-stl
import matplotlib.pyplot as plt # pip install matplotlib
import numpy as np # pip install numpy

# %%
script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
os.chdir(script_dir) # change the working directory

# load the STL file
geometry = mesh.Mesh.from_file('../sim/data/39_2-1-1-ReLA_212_thick_smooth.stl')

all_vertices = geometry.vectors.reshape(-1, 3) # flatten the triangle vertices into a (n*3, 3) array
all_vertices_unique = np.unique(all_vertices, axis=0)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(all_vertices_unique[:, 0], all_vertices_unique[:, 1], all_vertices_unique[:, 2], color='blue', s=1, marker='.')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_box_aspect([1,1,1])
plt.show()

print("min(xyz):")
print(np.min(all_vertices_unique, axis=0))
print("max(xyz):")
print(np.max(all_vertices_unique, axis=0))
print("mean(xyz):")
print(np.mean(all_vertices_unique, axis=0))

# %%