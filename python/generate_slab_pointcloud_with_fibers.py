"""
Create a sphere point cloud and compute fiber and sheet directions for each point.
This mirrors the C++ ComputeFiberAndSheetDirections math:
- face_norm is the surface normal (for a sphere it's radial: (pos - center)/|pos-center|)
- center_norm is the normalized position vector (same for centered sphere)
- circumferential direction = center_line x face_norm
- rotation angle beta = (beta_epi - beta_endo) * phi + beta_endo
- apply Rodrigues rotation of circumferential direction around face_norm by beta

Outputs a CSV with columns: x,y,z,phi,fiber_x,fiber_y,fiber_z,sheet_x,sheet_y,sheet_z

Usage:
    python3 generate_slab_pointcloud_with_fibers.py out.csv --nx 50 --ny 50 --nz 50 --dp 2.0 --radius 20.0

"""
#!/usr/bin/env python3
from types import SimpleNamespace
import math
import csv
import sys
from typing import Tuple
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Helper linear algebra utilities

def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-15:
        return np.zeros_like(v)
    return v / n


def rodrigues_rotate(v: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate vector v around unit axis by angle (radians) using Rodrigues' formula."""
    k = normalize(axis)
    if np.linalg.norm(k) < 1e-15:
        return np.zeros_like(v)
    return (math.cos(angle) * v + math.sin(angle) * np.cross(k, v) +
            (1.0 - math.cos(angle)) * np.dot(k, v) * k)


# Compute phi for a slab: normalized transmural coordinate (0..1)
# We'll treat slab thickness along y, with y in [ymin, ymax]
# phi = (y - ymin) / (ymax - ymin), but the C++ used something like probing diffusion field; this is an approximation

def compute_phi_from_slab_y(y: float, ymin: float, ymax: float) -> float:
    if ymax <= ymin:
        return 0.0
    return float((y - ymin) / (ymax - ymin))


def compute_phi_from_sphere_r(r: float, radius: float) -> float:
    # radial transmural coordinate from center (0) to surface (1)
    if radius <= 0.0:
        return 0.0
    return min(max(r / radius, 0.0), 1.0)


def compute_fiber_sheet_for_point(pos: np.ndarray,
                                 center_line_vector: np.ndarray,
                                 beta_epi: float,
                                 beta_endo: float,
                                 phi: float,
                                 center: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # For sphere: face_norm is the radial outward normal from the sphere center
    disp = pos - center
    r = np.linalg.norm(disp)
    if r < 1e-15:
        # at center, choose arbitrary normal
        face_norm = np.array([0.0, 1.0, 0.0])
    else:
        face_norm = disp / r

    center_norm = normalize(pos - center)

    # Ensure face_norm points outward relative to center_norm (mirror C++ flip logic)
    if np.dot(face_norm, center_norm) <= 0.0:
        face_norm = -face_norm

    circumferential_direction = np.cross(center_line_vector, face_norm)
    cd_norm = normalize(circumferential_direction)
    if np.linalg.norm(cd_norm) < 1e-15:
        # pick an arbitrary orthogonal vector
        if abs(face_norm[0]) < 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        else:
            tmp = np.array([0.0, 0.0, 1.0])
        cd_norm = normalize(np.cross(tmp, face_norm))

    beta = (beta_epi - beta_endo) * phi + beta_endo
    f_0 = rodrigues_rotate(cd_norm, face_norm, beta)
    f_0 = normalize(f_0)
    s_0 = face_norm.copy()
    return f_0, s_0

    # A better approach would be to set face_norm = [0, sign(y - center_y), 0]
    # But to match the C++ logic more generally, compute face_norm from displacement -> here assume face_norm = [0, +/-1, 0]

    # Decide face normal direction by pointing outward from slab center along y.
    # We'll choose slab center at y=0.0 by default: face_norm points from interior to exterior
    center_norm = normalize(pos)

    # For slab we take face_norm as unit vector along +y (epicardium outward direction).
    # But to mirror the C++ check (face_norm.dot(center_norm) <= 0 flip), we pick face_norm as +y first.
    face_norm = np.array([0.0, 1.0, 0.0])
    if np.dot(face_norm, center_norm) <= 0.0:
        face_norm = -face_norm

    # circumferential_direction = center_line_vector x face_norm
    circumferential_direction = np.cross(center_line_vector, face_norm)
    cd_norm = normalize(circumferential_direction)
    # If cd_norm is zero (centerline parallel to face_norm), pick an arbitrary orthogonal vector
    if np.linalg.norm(cd_norm) < 1e-15:
        # choose something orthogonal to face_norm
        if abs(face_norm[0]) < 0.9:
            tmp = np.array([1.0, 0.0, 0.0])
        else:
            tmp = np.array([0.0, 0.0, 1.0])
        cd_norm = normalize(np.cross(tmp, face_norm))

    # rotation angle
    beta = (beta_epi - beta_endo) * phi + beta_endo

    # Rodrigues' rotation: rotate cd_norm around face_norm by beta
    f_0 = rodrigues_rotate(cd_norm, face_norm, beta)

    # For slab we follow same interior cutoff: C++ sets fiber only for particles with pos.y < -ReferenceSpacing
    # Here, we won't zero fibers — the user can filter by y if needed. We normalize f_0.
    f_0 = normalize(f_0)
    s_0 = face_norm.copy()
    return f_0, s_0


# Explicit defaults so main() can be called without CLI arguments.
# Edit these values directly if you want different defaults.
args = SimpleNamespace(
    out_csv='out_sphere.csv',
    nx=40,
    ny=40,
    nz=40,
    dp=2.0,
    xmin=-35.0,
    xmax=35.0,
    ymin=-35.0,
    ymax=35.0,
    zmin=-35.0,
    zmax=35.0,
    beta_epi_deg=-70.0,
    beta_endo_deg=80.0,
    centerline=[0.0, 1.0, 0.0],
    radius=20.0,
    center=[0.0, 0.0, 0.0],
    plot=True,
    save_plot='',
    subsample=10,
)

beta_epi = math.radians(args.beta_epi_deg)
beta_endo = math.radians(args.beta_endo_deg)
center_line_vector = np.array(args.centerline, dtype=float)
center_line_vector = normalize(center_line_vector)

xs = np.linspace(args.xmin, args.xmax, args.nx)
ys = np.linspace(args.ymin, args.ymax, args.ny)
zs = np.linspace(args.zmin, args.zmax, args.nz)

rows = []
center = np.array(args.center, dtype=float)
for xi in xs:
    for yi in ys:
        for zi in zs:
            pos = np.array([xi, yi, zi], dtype=float)
            # include only points inside the sphere
            r = np.linalg.norm(pos - center)
            if r > args.radius:
                continue
            phi = compute_phi_from_sphere_r(r, args.radius)
            fiber, sheet = compute_fiber_sheet_for_point(pos, center_line_vector, beta_epi, beta_endo, phi, center)
            rows.append((xi, yi, zi, phi, fiber[0], fiber[1], fiber[2], sheet[0], sheet[1], sheet[2]))

# write CSV header and rows
with open(args.out_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['x', 'y', 'z', 'phi', 'fiber_x', 'fiber_y', 'fiber_z', 'sheet_x', 'sheet_y', 'sheet_z'])
    for r in rows:
        writer.writerow(["{:.9f}".format(v) if isinstance(v, float) else v for v in r])

print(f"Wrote {len(rows)} points to {args.out_csv}")

# Optional quick plotting
# Convert rows to arrays
data = np.array(rows)
if data.size == 0:
    print('No points to plot.')

X = data[:, 0].astype(float)
Y = data[:, 1].astype(float)
Z = data[:, 2].astype(float)
FX = data[:, 4].astype(float)
FY = data[:, 5].astype(float)
FZ = data[:, 6].astype(float)
PHI = data[:, 3].astype(float)

# subsample by stride to avoid overcrowding
stride = max(1, int(args.subsample))
ids = np.arange(0, X.size, stride)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
# scale arrows by dp for visibility
arrow_len = args.dp * 0.8
ax.quiver(X[ids], Y[ids], Z[ids], FX[ids], FY[ids], FZ[ids], length=arrow_len, normalize=True, cmap='viridis')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Fiber directions (subsampled)')
plt.show()

if args.save_plot:
    plt.savefig(args.save_plot, dpi=200)
    print(f'Saved plot to {args.save_plot}')
