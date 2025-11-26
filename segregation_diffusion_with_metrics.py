import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import os

# Reproducible seed (from your notes)
np.random.seed(0)

# Parameters (tweakable)
grid_size = 50
base_C = 0.5
sigma_C = 0.10
Zn_fraction = 0.05

D0 = 0.1               # base diffusion coefficient
Zn_D_factor = 0.25     # D at Zn sites = D0 * Zn_D_factor
h = 1.0                # grid spacing
dt = 0.1               # time-step (must satisfy dt <= h^2/(4*D0))
n_steps = 1000
report_every = 10

# Stability check (informational)
dt_max = h*h / (4.0 * D0)
print(f"Stability limit dt_max = {dt_max:.4g}; using dt = {dt}")

# Initialization
C = base_C + sigma_C * np.random.randn(grid_size, grid_size)
C = np.clip(C, 0.0, 1.0)
Z = (np.random.rand(grid_size, grid_size) < Zn_fraction).astype(int)
D = D0 * np.where(Z == 1, Zn_D_factor, 1.0)

# Utility: discrete 5-point Laplacian using periodic or Neumann boundary (here use Neumann via pad/reflection)
def laplacian(u):
    # Using numpy roll gives periodic BC. If you prefer Neumann (zero-flux) use mirrored padding or explicit handling.
    return (np.roll(u, -1, axis=0) + np.roll(u, 1, axis=0) +
            np.roll(u, -1, axis=1) + np.roll(u, 1, axis=1) - 4.0 * u) / (h * h)

# Moran's I (4-neighbor weights)
def morans_I(field):
    n = field.size
    x = field.ravel()
    xmean = x.mean()
    xd = x - xmean
    denom = np.sum(xd * xd)
    if denom == 0:
        return np.nan
    # build numerator with 4-neighbor adjacency via shifts
    W_sum = 0
    num = 0.0
    # neighbor shifts
    for shift in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        neighbor = np.roll(field, shift=shift, axis=(0,1)).ravel()
        num += np.sum(xd * (neighbor - xmean))
        W_sum += n  # each shift contributes n weights of 1
    I = (n / W_sum) * (num / denom)
    return I

# cluster stats: threshold-based (higher-than-threshold clusters)
def cluster_stats(field, threshold):
    mask = field > threshold
    labeled, ncomp = ndimage.label(mask, structure=np.array([[0,1,0],[1,1,1],[0,1,0]]))
    if ncomp == 0:
        return {"ncomp":0,"sizes":[], "max_size":0, "mean_size":0}
    sizes = ndimage.sum(mask, labeled, index=np.arange(1, ncomp+1))
    sizes = np.array(sizes, dtype=int)
    return {"ncomp": int(ncomp), "sizes": sizes, "max_size": int(sizes.max()), "mean_size": float(sizes.mean())}

# storage for metrics
metrics = {"step": [], "mean": [], "std": [], "cv": [], "moransI": [], "nclusters": [], "max_cluster_size": []}

# Simulation loop
for step in range(n_steps+1):
    if step % report_every == 0:
        meanC = C.mean()
        stdC = C.std(ddof=0)
        cv = stdC/meanC if meanC != 0 else np.nan
        I = morans_I(C)
        thr = meanC + 0.5 * stdC  # example threshold
        cs = cluster_stats(C, thr)
        metrics["step"].append(step)
        metrics["mean"].append(meanC)
        metrics["std"].append(stdC)
        metrics["cv"].append(cv)
        metrics["moransI"].append(I)
        metrics["nclusters"].append(cs["ncomp"])
        metrics["max_cluster_size"].append(cs["max_size"])
        print(f"step={step:5d} mean={meanC:.4f} std={stdC:.4f} CV={cv:.3f} MoranI={I:.4f} ncomp={cs['ncomp']} maxsize={cs['max_size']}")

    # Explicit diffusion update
    L = laplacian(C)
    C = C + (D * dt) * L
    # Optional: enforce bounds
    C = np.clip(C, 0.0, 1.0)

# Save outputs
outdir = "segregation_step_outputs"
os.makedirs(outdir, exist_ok=True)
np.savez_compressed(os.path.join(outdir, "step_arrays.npz"), C=C, Z=Z, D=D, metrics=metrics)

# Plots
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
im0 = axes[0,0].imshow(C, cmap='viridis')
axes[0,0].set_title("Final Cu concentration")
plt.colorbar(im0, ax=axes[0,0])
im1 = axes[0,1].imshow(Z, cmap='gray')
axes[0,1].set_title("Zn mask")
plt.colorbar(im1, ax=axes[0,1])

axes[1,0].plot(metrics["step"], metrics["std"], label="std")
axes[1,0].plot(metrics["step"], metrics["mean"], label="mean")
axes[1,0].set_xlabel("step")
axes[1,0].legend()
axes[1,0].set_title("Mean and Std vs time")

axes[1,1].plot(metrics["step"], metrics["moransI"], label="Moran's I")
axes[1,1].plot(metrics["step"], metrics["max_cluster_size"], label="max cluster size")
axes[1,1].legend()
axes[1,1].set_title("Spatial metrics vs time")

plt.tight_layout()
plt.savefig(os.path.join(outdir, "step_visuals.png"), dpi=150)
plt.show()
print("Saved outputs to", outdir)