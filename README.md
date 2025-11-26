# Sharma_Sagarika.
Development of Stable Cu-Based Argyrodites for Highly Efficient Energy Conversions
# Development of Stable Cu-Based Argyrodites for Highly Efficient Energy Conversions

**Researcher:** Sagarika Sharma  
**Project No.:** 10  
**Programming Method Used:** Python  

---



## 1. Background and Key Concepts

### 1.1 Overview

This report summarizes the current progress on developing a toy-model simulation to study Cu segregation and the stabilizing influence of Zn substitution in Cu-based argyrodites (Cu₆PS₅X). It covers milestones completed, the first step executed, challenges encountered, and next planned actions.

### 1.2 Problem: Cu Segregation
- During annealing (heat treatment), mobile Cu ions can **cluster**, leading to inhomogeneous regions.  
- Consequences:
  - Fluctuating electrical conductivity  
  - Unstable carrier concentration  
  - Poor reproducibility in experiments  

### 1.3 Zn Substitution as a Solution
- Introducing Zn²⁺ into some Cu sites:
  - Reduces Cu mobility locally  
  - Raises Cu vacancy formation energy  
  - Stabilizes the lattice  
- Expected effect: **lower segregation**, more uniform Cu distribution, more reproducible transport properties.

### 1.4 Key Metrics
- **Standard deviation (std)** of Cu concentration → primary segregation metric  
- Optional metrics: coefficient of variation, spatial autocorrelation (Moran’s I), cluster size

---

## 2. . Milestones Achieved
Milestone 0 – Setup

Status: Completed

0.1 Python Environment Installed

Python successfully installed (v3.x).

Required libraries (numpy, matplotlib) verifie.

### 2.1 Set reproducible randomness

np.random.seed(0) so results can be reproduced exactly.


import numpy as np
import matplotlib.pyplot as plt

# Step 1: Initialize grid with Cu concentration and Zn mask

# Reproducibility
np.random.seed(0)

### 2.2 Defined parameters

Grid size: 50×50 (2,500 cells).

Baseline Cu fraction base_C = 0.5.

Initial noise sigma_C = 0.10 (so each cell is 0.5 ± 0.1 sampled from a normal distribution).

Zn substitution target fraction: Zn_fraction = 0.05 (5% of cells chosen at random to contain Zn).

# Grid size
grid_size = 50

# Initial Cu concentration (random around 0.5)
C = 0.5 + 0.1 * np.random.randn(grid_size, grid_size)
C = np.clip(C, 0, 1)  # ensure values stay between 0 and 1

# Zn substitution mask (5%)
Zn_fraction = 0.05
Z = (np.random.rand(grid_size, grid_size) < Zn_fraction).astype(int)

### 2.3  Initialized Cu concentration field

C = base_C + sigma_C * np.random.randn(grid_size, grid_size) and clipped to [0,1].

### 2.4 Created Zn mask

Binary mask Z where each cell independently has probability 0.05 to be Zn.

Counted actual Zn sites and computed the actual fraction (since sampling gives an empirical fraction close to target).

# Plot initial Cu distribution and Zn mask
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

im1 = axes[0].imshow(C, cmap='viridis')
axes[0].set_title("Initial Cu Concentration")
plt.colorbar(im1, ax=axes[0])

im2 = axes[1].imshow(Z, cmap='viridis')
axes[1].set_title("Zn Mask (5%)")
plt.colorbar(im2, ax=axes[1])

plt.tight_layout()
plt.show()

## 3. Numeric results from this run

Grid size: 50 × 50 (2,500 cells)

Cu concentration:

Mean = 0.4981

Std = 0.0976

Min = 0.1883

Max = 0.8171

Zn mask:

Target fraction = 5.00%

Actual Zn count = 134 cells

Actual fraction = 5.36%

The actual Zn fraction is close to the target because of random sampling.

## 3.1 Saved outputs

Saved arrays C and Z to /mnt/data/segregation_step1_outputs/step1_arrays.npz.

Saved visualization PNG to /mnt/data/segregation_step1_outputs/step1_visuals.png.

## 3.2 Plotted results

Left: heatmap of initial Cu concentration (C).

Middle: Zn mask (binary).

Right: histogram of Cu fraction values across the grid.

## Files I saved


 
 








---

## 5. Interpretation

The Cu field is roughly centered at 0.5 with spread ≈ 0.1, as intended. The histogram confirms a near-Gaussian distribution clipped to [0,1].

Zn sites are sparse (~5%) and appear randomly scattered — this mimics a dilute substitution strategy in the toy model.

---

## 6. Suggested immediate next steps

Implement the diffusion update (the Laplacian-based smoothing) with spatially varying diffusion coefficient:

D_ij = D0 if no Zn, D0 * 0.25 if Zn present.

Use a standard 5-point Laplacian finite difference: L(C)[i,j] = C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1] - 4*C[i,j].

Choose D0 and a time-step dt such that the update C += D_ij * dt * Laplacian is stable (or use small dt).
