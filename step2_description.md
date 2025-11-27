# Step 2: Diffusion Simulation (Detailed Description)

## Overview

Step 2 simulates Cu ion diffusion in a 2D grid representing the Cu-based
argyrodite microstructure. The model applies a diffusion-like update
rule modified by the presence of Zn, which reduces local mobility.

## Goals

-   Demonstrate Cu segregation reduction with Zn substitution.
-   Track the standard deviation of Cu concentration over simulated
    annealing.
-   Save intermediate states for visualization and analysis.

## Model Components

### 1. Grid

-   Size: 50×50
-   `C[i,j]`: local Cu concentration (0--1)
-   `Z[i,j]`: Zn mask (0 or 1)
-   Initial Cu field generated with Gaussian noise.

### 2. Diffusion Equation

The update rule:

    C_{t+1} = C_t + D_{i,j} * Laplacian(C_t) * dt

with periodic boundary conditions.

### 3. Diffusion Coefficient

    D[i,j] = D0             if Z[i,j] == 0
    D[i,j] = D0 * 0.25      if Z[i,j] == 1

### 4. Simulation Parameters

-   `D0 = 0.2`
-   Zn fractions tested: **0%, 5%, 10%**
-   `dt = 0.2`
-   Total steps: **150**
-   Snapshots saved at steps 0, 50, 100, 150

### 5. Outputs

-   Final Cu distribution `C_final`
-   Time series `std(C)`
-   `.npz` data files for each Zn fraction
-   Plot: **std_vs_time.png**

## Interpretation

-   Zn reduces Cu mobility locally.
-   Zn → faster homogenization → reduced segregation.
-   Standard deviation of Cu decreases more rapidly with Zn.

## Files Produced

-   `step2_z0.npz`
-   `step2_z5.npz`
-   `step2_z10.npz`
-   `std_vs_time.png`
