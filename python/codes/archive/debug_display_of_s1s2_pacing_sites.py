import numpy as np
import plotly.graph_objects as go # pip install plotly. For 3D interactive plot: triangular mesh, and activation movie

def execute(voxel, s1_pacing_voxel_id, s2_pacing_voxel_id):
    # Extract voxel groups
    voxel_s1 = voxel[s1_pacing_voxel_id, :]
    voxel_s2 = voxel[s2_pacing_voxel_id, :]

    id = np.arange(voxel.shape[0])
    id = id[~np.isin(id, s1_pacing_voxel_id)]
    id = id[~np.isin(id, s2_pacing_voxel_id)]
    voxel_rest = voxel[id, :]

    # Create 3D scatter plot
    fig = go.Figure()

    # S1 pacing voxels (blue)
    fig.add_trace(go.Scatter3d(
        x=voxel_s1[:,0], y=voxel_s1[:,1], z=voxel_s1[:,2],
        mode='markers',
        marker=dict(size=3, color='blue'),
        name='S1 pacing'
    ))

    # S2 pacing voxels (red)
    fig.add_trace(go.Scatter3d(
        x=voxel_s2[:,0], y=voxel_s2[:,1], z=voxel_s2[:,2],
        mode='markers',
        marker=dict(size=3, color='red'),
        name='S2 pacing'
    ))

    # Remaining voxels (grey)
    fig.add_trace(go.Scatter3d(
        x=voxel_rest[:,0], y=voxel_rest[:,1], z=voxel_rest[:,2],
        mode='markers',
        marker=dict(size=2, color='grey'),
        name='Other voxels'
    ))

    # Layout
    fig.update_layout(
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z"
        ),
        width=800,
        height=700,
        title="Voxel 3D Plot"
    )

    fig.show()
