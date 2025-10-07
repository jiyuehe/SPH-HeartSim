#%%
import os, sys
import codes
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

#%%
class MeshSelector:
    def __init__(self, vertices, faces, data_path, vertex_flag, vertex_color):
        self.vertices = vertices
        self.faces = faces
        self.data_path = data_path
        self.vertex_flag = vertex_flag
        self.vertex_color = vertex_color
        
        # Define color map for different flags
        color_map = ['red', 'green', 'blue', 'yellow', 'cyan', 'magenta', 'orange', 'purple']
        self.color_map = color_map

        self.selection_polygon = []
        self.is_selecting = False
        self.selection_mode = False # Toggle between rotate and select modes
        
        # Create figure and 3D axis
        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        # Plot the mesh
        self.mesh = self.ax.plot_trisurf(
            self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2],
            triangles=self.faces,
            linewidth=0.2, edgecolor='gray', alpha=0, color='white'
        )

        self.flagged_scatter = self.ax.scatter([], [], [], c=[], s=5, depthshade=False, marker='.')

        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        codes.set_axes_equal.execute(self.ax)
        self.ax.set_title('3D Mesh Vertex Selection Tool')
        self.update_display()

        # Connect mouse events
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        
        # Add text input box for flag value
        from matplotlib.widgets import Button, TextBox
        ax_textbox = plt.axes([0.15, 0.90, 0.1, 0.04])
        self.textbox = TextBox(ax_textbox, 'Flag Value:', initial='1')
        
        # Add save button
        ax_save = plt.axes([0.81, 0.90, 0.08, 0.04])
        self.btn_save = Button(ax_save, 'Save')
        self.btn_save.on_clicked(self.save_selection)

        # Instructions
        self.fig.text(0.5, 0.02, 
            'Press "a" to toggle Rotate/Select mode | Mouse left drag: Add vertices to selection | Right click: Clear all selected',
            ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        self.mode_text = self.fig.text(0.5, 0.95, 'Mode: Rotate (press "a" change to Select)', 
            ha='center', fontsize=12, weight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        plt.show()

    def update_display(self):
        # Highlight selected vertices (where flag > 0) with different colors
        flagged_idx = np.where(self.vertex_flag > 0)[0]
        if len(flagged_idx) == 0:
            self.flagged_scatter._offsets3d = ([], [], [])
            self.flagged_scatter.set_facecolor([])
        else:
            verts = self.vertices[flagged_idx]
            print(verts)
            colors = []
            for idx in flagged_idx:
                colors.append(self.color_map[(self.vertex_flag[idx]-1) % len(self.color_map)])
            self.flagged_scatter._offsets3d = (verts[:,0], verts[:,1], verts[:,2])
            self.flagged_scatter.set_color(colors)

    def on_press(self, event):
        # Handle mouse press events
        if event.inaxes != self.ax:
            return
        
        if event.button == 1 and self.selection_mode: # Left click in selection mode
            self.is_selecting = True
            self.selection_polygon = [(event.xdata, event.ydata)]
        elif event.button == 3: # Right click - clear all selections
            self.vertex_flag = np.zeros(len(self.vertices), dtype=int)
            self.selection_polygon = []
            self.update_display()
            print("All selections cleared")
    
    def on_motion(self, event):
        # If you’re not currently dragging (i.e., is_selecting is False) 
        # or the mouse isn’t inside the 3D axes, do nothing.
        if not self.is_selecting or event.inaxes != self.ax:
            return
        
        # Record the mouse position
        if event.xdata is not None and event.ydata is not None:
            self.selection_polygon.append((event.xdata, event.ydata))
        
    def on_release(self, event):
        # Handle mouse release - complete selection
        if not self.is_selecting or event.button != 1:
            return
        
        self.is_selecting = False
        
        if len(self.selection_polygon) < 3:
            print("Selection too small, need at least 3 points")
            self.selection_polygon = []
            return
        
        # Close the polygon
        self.selection_polygon.append(self.selection_polygon[0])
        
        # Project vertices to 2D screen coordinates and update flags
        self.select_vertices_in_polygon()
        
        # Clear polygon and redraw
        self.selection_polygon = []
        self.update_display()
    
    def on_key(self, event):
        # Handle key press events
        if event.key.lower() == 'a':
            self.selection_mode = not self.selection_mode
            if self.selection_mode:
                self.mode_text.set_text('Mode: Select (press "a" change to Rotate)')
                self.mode_text.set_bbox(dict(boxstyle='round', facecolor='yellow', alpha=0.8))
                # Disable 3D rotation
                self.ax.disable_mouse_rotation()
                print("Selection mode ON - Click and drag to add vertices to selection")
            else:
                self.mode_text.set_text('Mode: Rotate (press "a" change to Select)')
                self.mode_text.set_bbox(dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
                # Enable 3D rotation
                self.ax.mouse_init()
                print("Rotation mode ON - Drag to rotate view")
            self.fig.canvas.draw_idle()
    
    def save_selection(self, event):
        # Save vertex flag array to a numpy file
        num_selected = np.sum(self.vertex_flag)
        if num_selected == 0:
            print("No vertices selected to save!")
            return
        
        filename = 'vertex_flag.npy'
        np.save(self.data_path + filename, self.vertex_flag)
        print(f"Saved vertex flags to '{filename}' ({num_selected} vertices selected)")
    
    def select_vertices_in_polygon(self):
        # Select vertices that fall within the drawn polygon (front face only)
        # and set their flags to the current flag value (accumulative)
        if len(self.selection_polygon) < 3:
            return
        
        # Get the flag value from text box
        flag_value = int(self.textbox.text)
        
        # Get the current 3D to 2D projection matrix
        proj_matrix = self.ax.get_proj()
        
        # Transform 3D vertices to 2D display coordinates and get depths
        vertices_2d = []
        depths = []
        for vertex in self.vertices:
            vec = np.array([vertex[0], vertex[1], vertex[2], 1.0])
            proj = proj_matrix @ vec
            
            if proj[3] != 0:
                x_2d = proj[0] / proj[3]
                y_2d = proj[1] / proj[3]
                z_depth = proj[2] / proj[3]
            else:
                x_2d = proj[0]
                y_2d = proj[1]
                z_depth = proj[2]
            
            vertices_2d.append([x_2d, y_2d])
            depths.append(z_depth)
        
        vertices_2d = np.array(vertices_2d)
        depths = np.array(depths)
        
        # Create a path from the selection polygon
        path = Path(self.selection_polygon)
        
        # Check which vertices are inside the polygon
        inside = path.contains_points(vertices_2d)
        
        # For each vertex inside, check if it's visible (front-facing)
        inside_indices = np.where(inside)[0]
        if len(inside_indices) > 0:
            min_depth = np.min(depths[inside_indices])
            tol = 1e-3  # tolerance for front-face detection
            front_facing = np.abs(depths - min_depth) < tol
            inside = inside & front_facing
        
        # Update vertex flags for newly selected vertices
        newly_selected = np.sum(inside & (self.vertex_flag == 0))
        self.vertex_flag[inside] = flag_value
        
        total_selected = np.sum(self.vertex_flag > 0)
        print(f"Set {np.sum(inside)} vertices to flag {flag_value} (total selected: {total_selected})")

#%%
if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
    os.chdir(script_dir) # change the working directory

    folder_path = "../build/sim/bin/output/"
    t, voltage, gate_variable, stress, xyz = codes.load_simulation_result.execute(folder_path)
    # t[time_id]
    # voltage[particles, time_steps]
    # stress[particles, time_steps]
    # xyz[particles, time_steps, coordinates]

    particle = xyz[:,0,:] # particle coordinates at time 0

    if os.path.exists(folder_path + 'particle_flag.npy'): # file exist
        particle_flag = np.load(folder_path + 'particle_flag.npy')
    else: # file do not exist
        particle_flag = np.zeros(len(particle), dtype=int)

    particle_color = np.ones((len(particle_flag), 3))  # shape (11, 3), all ones

    selector = MeshSelector(vertex, face, data_path, vertex_flag, vertex_color)

# %%
