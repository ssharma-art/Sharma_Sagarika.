"""
cu_diffusion_sim.py

Simple Cu diffusion toy model (NumPy).
Compares "No Zn" vs "With Zn" by reducing local diffusion at Zn sites.

Requirements:
    pip install numpy matplotlib

Run:
    python cu_diffusion_sim.py
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Parameters (easy to change)
# ----------------------------
nx, ny = 60, 60           # grid size
steps = 300               # time steps to simulate
save_steps = [0, steps-1] # snapshots to show (start and end)
zn_fraction = 0.05        # fraction of sites that are Zn (5%)
D_base = 0.18             # base diffusion factor
zn_D_factor = 0.25        # multiplier for D where Zn present (reduces mobility)
noise_amp = 0.02          # initial random noise amplitude
seed = 2025               # random seed for reproducibility

np.random.seed(seed)

# ----------------------------
# Initialization
# ----------------------------
# Base Cu concentration (mean) and small random fluctuations
Cu0 = 0.6 + noise_amp * (np.random.rand(nx, ny) - 0.5)

# Zn masks for two cases
Zn_none = np.zeros((nx, ny), dtype=bool)
Zn_some = (np.random.rand(nx, ny) < zn_fraction)

# Two independent simulations: copy initial Cu
Cu_nozn = Cu0.copy()
Cu_wzn = Cu0.copy()

# Arrays to store std over time
std_nozn = np.zeros(steps)
std_wzn = np.zeros(steps)

# ----------------------------
# Diffusion update function
# ----------------------------
def diffusion_step(Cu, Zn_mask, D_base=D_base, zn_factor=zn_D_factor):
    """
    One explicit diffusion-like update using 4-neighbor Laplacian.
    Cu : 2D array
    Zn_mask : boolean 2D array, same shape
    Returns new Cu array (clipped to [0,1])
    """
    # periodic boundary neighbors (roll wraps around)
    up = np.roll(Cu, 1, axis=0)
    down = np.roll(Cu, -1, axis=0)
    left = np.roll(Cu, 1, axis=1)
    right = np.roll(Cu, -1, axis=1)
    lap = (up + down + left + right - 4.0 * Cu)

    # local diffusion coefficient (reduced at Zn sites)
    D_local = D_base * (np.where(Zn_mask, zn_factor, 1.0))

    Cu_new = Cu + D_local * lap
    # keep physically bounded
    return np.clip(Cu_new, 0.0, 1.0)

# ----------------------------
# Run simulation
# ----------------------------
for t in range(steps):
    Cu_nozn = diffusion_step(Cu_nozn, Zn_none)
    Cu_wzn = diffusion_step(Cu_wzn, Zn_some)

    std_nozn[t] = np.std(Cu_nozn)
    std_wzn[t] = np.std(Cu_wzn)

# ----------------------------
# Plots: std vs time
# ----------------------------
plt.figure(figsize=(7,4))
plt.plot(std_nozn, label="No Zn")
plt.plot(std_wzn, label=f"With Zn ({int(zn_fraction*100)}%)")
plt.xlabel("Step")
plt.ylabel("Std(Cu) — segregation metric")
plt.title("Cu segregation over time (toy diffusion model)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("std_vs_time.png", dpi=200)
print("Saved std_vs_time.png")

# ----------------------------
# Plots: Cu maps (initial and final)
# ----------------------------
fig, axes = plt.subplots(2, 2, figsize=(8,7))

# initial Cu map (same for both sims)
axes[0,0].imshow(Cu0, vmin=0.45, vmax=0.75)
axes[0,0].set_title("Initial Cu (both cases)")
axes[0,0].axis('off')

axes[0,1].axis('off')  # empty cell for layout symmetry

# final maps
axes[1,0].imshow(Cu_nozn, vmin=0.45, vmax=0.75)
axes[1,0].set_title("Final Cu - No Zn")
axes[1,0].axis('off')

axes[1,1].imshow(Cu_wzn, vmin=0.45, vmax=0.75)
axes[1,1].set_title(f"Final Cu - With Zn ({int(zn_fraction*100)}%)")
axes[1,1].axis('off')

plt.suptitle("Cu concentration maps (toy model)")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("cu_maps.png", dpi=200)
print("Saved cu_maps.png")
plt.show()

# ----------------------------
# Final numeric summary
# ----------------------------
print("Initial std(Cu): {:.6f}".format(np.std(Cu0)))
print("Final std(Cu) - No Zn: {:.6f}".format(std_nozn[-1]))
print("Final std(Cu) - With Zn: {:.6f}".format(std_wzn[-1]))
