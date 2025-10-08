#%%
import os, sys
import codes
import numpy as np # pip install numpy
import matplotlib.pyplot as plt # pip install matplotlib
from matplotlib.path import Path

#%%
class MeshSelector:
    def __init__(self, node, node_flag, folder_path):
        self.node = node
        self.node_flag = node_flag
        self.folder_path = folder_path
        
        # Define color map for different flags
        color_map = ['red', 'green', 'blue', 'yellow', 'cyan', 'magenta', 'orange', 'purple']
        self.color_map = color_map

        self.selection_polygon = []
        self.is_selecting = False
        self.selection_mode = False # Toggle between rotate and select modes
        
        # Create figure and 3D axis
        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')

        self.scatter = self.ax.scatter(self.node[:, 0], self.node[:, 1], self.node[:, 2], c='grey', s=5, depthshade=False, marker='.')

        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        codes.set_axes_equal.execute(self.ax)
        self.ax.set_title('3D Node Selection Tool')
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
            'Press "a" to toggle Rotate/Select mode | Mouse left drag: Add nodes to selection | Right click: Clear all selected',
            ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        self.mode_text = self.fig.text(0.5, 0.95, 'Mode: Rotate (press "a" change to Select)', 
            ha='center', fontsize=12, weight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        plt.show()

    def update_display(self):
        N = self.node.shape[0] # number of nodes
        colors = ['grey'] * N # color for those not selected

        # assign colors to different flags
        flagged_id = np.where(self.node_flag > 0)[0]
        for id in flagged_id:
            colors[id] = self.color_map[(self.node_flag[id] - 1) % len(self.color_map)]

        self.scatter.set_color(colors)

    def on_press(self, event):
        # Handle mouse press events
        if event.inaxes != self.ax:
            return
        
        if event.button == 1 and self.selection_mode: # Left click in selection mode
            self.is_selecting = True
            self.selection_polygon = [(event.xdata, event.ydata)]
        elif event.button == 3: # Right click - clear all selections
            self.node_flag = np.zeros(len(self.node), dtype=int)
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
        
        # Project node to 2D screen coordinates and update flags
        self.select_node_in_polygon()
        
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
                print("Selection mode ON - Click and drag to add nodes to selection")
            else:
                self.mode_text.set_text('Mode: Rotate (press "a" change to Select)')
                self.mode_text.set_bbox(dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
                # Enable 3D rotation
                self.ax.mouse_init()
                print("Rotation mode ON - Drag to rotate view")
            self.fig.canvas.draw_idle()
    
    def save_selection(self, event):
        # Save vertex flag array to a numpy file
        num_selected = np.sum(self.node_flag)
        if num_selected == 0:
            print("No node selected to save!")
            return
        
        filename = 'node_flag.npy'
        np.save(self.folder_path + filename, self.node_flag)
        print(f"Saved vertex flags to '{filename}' ({num_selected} nodes selected)")
    
    def select_node_in_polygon(self):
        # Select node that fall within the drawn polygon (front face only)
        # and set their flags to the current flag value (accumulative)
        if len(self.selection_polygon) < 3:
            return
        
        # Get the flag value from text box
        flag_value = int(self.textbox.text)
        
        # Get the current 3D to 2D projection matrix
        proj_matrix = self.ax.get_proj()
        
        # Transform 3D node to 2D display coordinates and get depths
        node_2d = []
        depths = []
        for vertex in self.node:
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
            
            node_2d.append([x_2d, y_2d])
            depths.append(z_depth)
        
        node_2d = np.array(node_2d)
        depths = np.array(depths)
        
        # Create a path from the selection polygon
        path = Path(self.selection_polygon)
        
        # Check which node are inside the polygon
        inside = path.contains_points(node_2d)
        
        # For each vertex inside, check if it's visible (front-facing)
        inside_indices = np.where(inside)[0]
        if len(inside_indices) > 0:
            min_depth = np.min(depths[inside_indices])
            tol = 1e-3  # tolerance for front-face detection
            front_facing = np.abs(depths - min_depth) < tol
            inside = inside & front_facing
        
        # Update vertex flags for newly selected node
        newly_selected = np.sum(inside & (self.node_flag == 0))
        self.node_flag[inside] = flag_value
        
        total_selected = np.sum(self.node_flag > 0)
        print(f"Set {np.sum(inside)} nodes to flag {flag_value} (total selected: {total_selected})")

#%%
if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__)) # get the path of the current script
    os.chdir(script_dir) # change the working directory

    folder_path = "../build/sim/bin/output/"
    t, voltage, gate_variable, stress, xyz = codes.load_simulation_result.execute(folder_path)
    # t[time_id]
    # voltage[nodes, time_steps]
    # stress[nodes, time_steps]
    # xyz[nodes, time_steps, coordinates]

    node = xyz[:,0,:] # node coordinates at time 0

    if os.path.exists(folder_path + 'node_flag.npy'): # file exist
        node_flag = np.load(folder_path + 'node_flag.npy')
    else: # file do not exist
        node_flag = np.zeros(len(node), dtype=int)

    np.set_printoptions(threshold=np.inf) # disable summarization, print out all elements

    s1_pacing_node_id = np.where(node_flag == 1)[0]
    print('s1 pacing nodes: ')
    print(", ".join(map(str, s1_pacing_node_id))) # "," in between elements

    s2_pacing_node_id = np.where(node_flag == 2)[0]
    print('s2 pacing nodes: ')
    print(", ".join(map(str, s2_pacing_node_id))) # "," in between elements

    selector = MeshSelector(node, node_flag, folder_path)

# %%
