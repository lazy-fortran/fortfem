#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.tri as tri

def load_data(filename):
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1], data[:, 2]

# Load data
x_lap, y_lap, z_lap = load_data("lapack_solution.dat")
x_umf, y_umf, z_umf = load_data("umfpack_solution.dat")
x_diff, y_diff, z_diff = load_data("difference.dat")

# Create triangulation
triang = tri.Triangulation(x_lap, y_lap)

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("2D Poisson Solver Comparison: LAPACK vs UMFPACK", fontsize=14)

# LAPACK solution
cs1 = axes[0,0].tricontourf(triang, z_lap, levels=20, cmap="viridis")
axes[0,0].set_title("LAPACK Solution")
axes[0,0].set_xlabel("x")
axes[0,0].set_ylabel("y")
axes[0,0].set_aspect("equal")
plt.colorbar(cs1, ax=axes[0,0])

# UMFPACK solution
cs2 = axes[0,1].tricontourf(triang, z_umf, levels=20, cmap="viridis")
axes[0,1].set_title("UMFPACK Solution")
axes[0,1].set_xlabel("x")
axes[0,1].set_ylabel("y")
axes[0,1].set_aspect("equal")
plt.colorbar(cs2, ax=axes[0,1])

# Difference
cs3 = axes[1,0].tricontourf(triang, z_diff, levels=20, cmap="plasma")
axes[1,0].set_title("Absolute Difference")
axes[1,0].set_xlabel("x")
axes[1,0].set_ylabel("y")
axes[1,0].set_aspect("equal")
plt.colorbar(cs3, ax=axes[1,0])

# Exact solution for comparison
z_exact = np.sin(np.pi * x_lap) * np.sin(np.pi * y_lap)
cs4 = axes[1,1].tricontourf(triang, z_exact, levels=20, cmap="viridis")
axes[1,1].set_title("Exact Solution")
axes[1,1].set_xlabel("x")
axes[1,1].set_ylabel("y")
axes[1,1].set_aspect("equal")
plt.colorbar(cs4, ax=axes[1,1])

plt.tight_layout()
plt.savefig("solver_comparison_2d.png", dpi=300, bbox_inches="tight")
plt.show()

# Print statistics
print(f"Max difference: {np.max(z_diff):.2e}")
print(f"RMS difference: {np.sqrt(np.mean(z_diff**2)):.2e}")
print(f"Max LAPACK solution: {np.max(np.abs(z_lap)):.6f}")
print(f"Max UMFPACK solution: {np.max(np.abs(z_umf)):.6f}")
print(f"Max exact solution: {np.max(np.abs(z_exact)):.6f}")
