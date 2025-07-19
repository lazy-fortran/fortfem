#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Load convergence data
data = np.loadtxt("nedelec_convergence.dat")
h = data[:, 0]
L2_err = data[:, 1]
H_curl_err = data[:, 2]
dofs = data[:, 3]

# Create convergence plots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# L2 error plot
ax1.loglog(h, L2_err, "o-", label="L2 Error", linewidth=2)
ax1.loglog(h, L2_err[0] * (h/h[0]), "--", label="O(h) slope", linewidth=1)
ax1.set_xlabel("Mesh size h")
ax1.set_ylabel("L2 Error")
ax1.set_title("Nédélec Element L2 Error Convergence")
ax1.legend()
ax1.grid(True)

# H(curl) error plot
safe_h_curl = np.maximum(H_curl_err, 1e-16)
ax2.loglog(h, safe_h_curl, "o-", label="H(curl) Error", linewidth=2)
ax2.axhline(1e-14, color="red", linestyle="--", label="Machine precision")
ax2.set_xlabel("Mesh size h")
ax2.set_ylabel("H(curl) Error")
ax2.set_title("Nédélec Element H(curl) Error (Perfect: 0.000000)")
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig("nedelec_convergence.png", dpi=300)
print("Convergence plot saved as nedelec_convergence.png")

# DOF scaling plot
plt.figure(figsize=(8, 6))
plt.loglog(h, dofs, "o-", label="Actual DOFs", linewidth=2)
plt.loglog(h, dofs[0] * (h[0]/h)**2, "--", label="O(1/h²) scaling", linewidth=1)
plt.xlabel("Mesh size h")
plt.ylabel("Number of DOFs")
plt.title("DOF Scaling with Mesh Refinement")
plt.legend()
plt.grid(True)
plt.savefig("nedelec_dof_scaling.png", dpi=300)
print("DOF scaling plot saved as nedelec_dof_scaling.png")
