

Step 3: Prototype diffusion simulation for Cu redistribution with Zn substitution.
- Implements explicit finite-difference diffusion with spatially-varying D.
- Computes metrics: std(C), coefficient of variation (CV), optional Moran's I.
- Saves arrays and figures for analysis.

Usage:
    python step3_prototype.py
or run the cells in a Jupyter notebook.

Author: (your name)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import json
from datetime import datetime

# -------------------------
# PARAMETERS (tweak these)
# -------------------------
params = {
    "grid_size": 50,         # L x L grid
    "base_C": 0.5,           # baseline Cu fraction
    "sigma_C": 0.10,         # initial perturbation std
    "Zn_fractions": [0.0, 0.05, 0.10],  # Zn substitution fractions to test
    "D0": 0.2,               # base diffusion coefficient
    "Zn_factor": 0.25,       # multiplier for D where Zn present (so D = D0 * Zn_factor)
    "dt": 0.2,               # time step
    "nsteps": 300,           # number of time steps
    "save_every": 50,        # save snapshots every N steps
    "ensemble_runs": 1,      # number of independent RNG seeds per Zn fraction (set >1 for statistics)
    "out_dir": "step3_outputs",
    "compute_moran": True,   # compute Moran's I (spatial autocorrelation) each step (costly)
}

# Create outputs folder
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
out_dir = os.path.join(params["out_dir"], timestamp)
os.makedirs(out_dir, exist_ok=True)

# Save params for provenance
with open(os.path.join(out_dir, "params.json"), "w") as f:
    json.dump(params, f, indent=2)

# -------------------------
# Utility functions
# -------------------------
def init_fields(grid_size, base_C, sigma_C, seed=None):
    """Initialize Cu concentration field C and return it (clipped to [0,1])."""
    rng = np.random.RandomState(seed)
    C0 = base_C + sigma_C * rng.randn(grid_size, grid_size)
    C0 = np.clip(C0, 0.0, 1.0)
    return C0, rng

def make_Z_mask(grid_size, zn_fraction, rng):
    """Return binary Zn mask Z with empirical fraction near zn_fraction."""
    Z = (rng.rand(grid_size, grid_size) < zn_fraction).astype(int)
    return Z

def laplacian_periodic(A):
    """5-point Laplacian with periodic boundary conditions using np.roll."""
    return (np.roll(A, 1, axis=0) + np.roll(A, -1, axis=0) +
            np.roll(A, 1, axis=1) + np.roll(A, -1, axis=1) - 4*A)

def compute_metrics(C):
    """Return (mean, std, cv)."""
    mean = C.mean()
    std = C.std(ddof=0)
    cv = std / mean if mean != 0 else np.nan
    return mean, std, cv

def morans_I(C):
    """
    Compute a simple global Moran's I for the 2D grid using rook adjacency (4 neighbors).
    This is a simple version (no row-standardization). For larger grids or many steps,
    this will cost more CPU.
    """
    n = C.size
    mean_C = C.mean()
    W = 0
    num = 0.0
    den = np.sum((C - mean_C)**2)
    # weight neighbors with 1 for 4 neighbors
    # use roll to iterate neighbors
    for shift in [(1,0),(-1,0),(0,1),(0,-1)]:
        neighbor = np.roll(np.roll(C, shift[0], axis=0), shift[1], axis=1)
        num += np.sum((C - mean_C) * (neighbor - mean_C))
        W += n  # each neighbor relation contributes n weights of 1
    if den == 0 or W == 0:
        return np.nan
    I = (n / W) * (num / den)
    return I

# -------------------------
# Single simulation runner
# -------------------------
def run_simulation(C0, Z, D0, Zn_factor, dt, nsteps, save_every, compute_moran=False):
    """Run simulation and return results dictionary."""
    grid_size = C0.shape[0]
    # spatially varying diffusion coefficient
    D = D0 * np.ones_like(C0)
    D[Z == 1] = D0 * Zn_factor

    C = C0.copy()
    times = []
    stds = []
    means = []
    cvs = []
    morans = []
    snapshots = {0: C.copy()}

    # stability check hint (dx = 1 assumed)
    maxD = D.max()
    stability_param = 4 * maxD * dt
    if stability_param > 1.0:
        print(f"WARNING: Explicit scheme may be unstable: 4*D_max*dt = {stability_param:.3f} > 1")
        # you may wish to reduce dt or D0

    for step in range(1, nsteps + 1):
        L = laplacian_periodic(C)
        C = C + dt * D * L
        C = np.clip(C, 0.0, 1.0)

        mean, std, cv = compute_metrics(C)
        times.append(step * dt)
        means.append(mean)
        stds.append(std)
        cvs.append(cv)
        if compute_moran:
            morans.append(morans_I(C))

        if step % save_every == 0:
            snapshots[step] = C.copy()

    results = {
        "C_final": C.copy(),
        "times": np.array(times),
        "means": np.array(means),
        "stds": np.array(stds),
        "cvs": np.array(cvs),
        "morans": np.array(morans) if compute_moran else None,
        "snapshots": snapshots,
        "D": D
    }
    return results

# -------------------------
# Top-level experiment loop (Zn fractions and ensembles)
# -------------------------
all_results = {}  # keyed by zn_fraction -> list of ensemble results
grid_size = params["grid_size"]

for zf in params["Zn_fractions"]:
    all_results[zf] = []
    for run_idx in range(params["ensemble_runs"]):
        seed = run_idx  # or use different seed scheme
        C0, rng = init_fields(grid_size, params["base_C"], params["sigma_C"], seed=seed)
        Z = make_Z_mask(grid_size, zf, rng)
        res = run_simulation(C0, Z, params["D0"], params["Zn_factor"],
                             params["dt"], params["nsteps"], params["save_every"],
                             compute_moran=params["compute_moran"])
        # add metadata
        res["C0"] = C0
        res["Z"] = Z
        res["seed"] = seed
        all_results[zf].append(res)

        # save numeric results per run
        fname = os.path.join(out_dir, f"results_z{int(zf*100)}_run{run_idx}.npz")
        np.savez_compressed(fname,
                            C0=C0, Z=Z, C_final=res["C_final"],
                            times=res["times"], means=res["means"], stds=res["stds"],
                            cvs=res["cvs"], morans=res["morans"] if res["morans"] is not None else np.array([]))
        print(f"Saved {fname} (Zn {zf}, run {run_idx})")

# -------------------------
# Plotting utilities (aggregate / individual)
# -------------------------
def plot_std_vs_time(all_results, out_path):
    plt.figure(figsize=(7,4))
    for zf, runs in all_results.items():
        # average across ensemble if multiple
        if len(runs) == 1:
            times = runs[0]["times"]
            stds = runs[0]["stds"]
            plt.plot(times, stds, label=f"Zn {zf*100:.0f}%")
        else:
            # compute ensemble mean and shading
            stack = np.stack([r["stds"] for r in runs], axis=0)
            mean = stack.mean(axis=0)
            std_err = stack.std(axis=0) / np.sqrt(len(runs))
            plt.plot(runs[0]["times"], mean, label=f"Zn {zf*100:.0f}% mean")
            plt.fill_between(runs[0]["times"], mean - std_err, mean + std_err, alpha=0.2)
    plt.xlabel("Time (arb. units)")
    plt.ylabel("Std(C)")
    plt.title("Std(C) vs time")
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    print("Saved:", out_path)

def plot_snapshot(C, title, path):
    plt.figure(figsize=(4,4))
    plt.imshow(C, origin='lower', vmin=0.0, vmax=1.0)
    plt.title(title)
    plt.colorbar(label="Cu fraction")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print("Saved:", path)

# Save std vs time plot
plot_std_vs_time(all_results, os.path.join(out_dir, "std_vs_time.png"))

# Save initial and final snapshots for each Zn fraction (first ensemble run)
for zf, runs in all_results.items():
    r = runs[0]
    plot_snapshot(r["C0"], f"Initial C (Zn {zf*100:.0f}%)", os.path.join(out_dir, f"init_z{int(zf*100)}.png"))
    plot_snapshot(r["C_final"], f"Final C (Zn {zf*100:.0f}%)", os.path.join(out_dir, f"final_z{int(zf*100)}.png"))
    # optionally save D map and Z mask
    plot_snapshot(r["D"], f"D map (Zn {zf*100:.0f}%)", os.path.join(out_dir, f"Dmap_z{int(zf*100)}.png"))
    plot_snapshot(r["Z"], f"Zn mask (Zn {zf*100:.0f}%)", os.path.join(out_dir, f"Zmask_z{int(zf*100)}.png"))

# -------------------------
# Save a summary JSON with key numbers for quick inspection
# -------------------------
summary = {}
for zf, runs in all_results.items():
    vals = []
    for r in runs:
        initial_std = r["C0"].std()
        final_std = r["C_final"].std()
        vals.append({"seed": r["seed"], "initial_std": float(initial_std), "final_std": float(final_std),
                     "zn_count": int(r["Z"].sum()), "zn_fraction_emp": float(r["Z"].mean())})
    summary[float(zf)] = vals

with open(os.path.join(out_dir, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("Summary saved to", os.path.join(out_dir, "summary.json"))

print("All outputs saved to:", out_dir)
