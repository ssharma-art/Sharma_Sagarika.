# Step 3: Python Prototype Implementation --- Detailed Description

## Overview

Step 3 is the core implementation stage of the simulation framework.\
Here, you translate the conceptual model (from Steps 1 and 2) into a
working Python prototype that performs:

1.  **Cu diffusion simulation with Zn-modified mobility**
2.  **Metric tracking (standard deviation and optional coefficient of
    variation)**
3.  **Data saving**
4.  **Visualization**

This step transforms the theoretical framework into an actionable
computational tool.

------------------------------------------------------------------------

# 1. Objectives of Step 3

The goals of this step are:

-   Implement the diffusion update rule numerically.
-   Incorporate Zn masks into the diffusion coefficient.
-   Run the simulation loop over many time steps.
-   Collect time-series data for later analysis.
-   Output data in `.npz` and figures in `.png`.
-   Validate visually that the simulation behaves as expected.

This establishes the foundation for all later scientific analysis.

------------------------------------------------------------------------

# 2. Code Structure

The prototype is divided into four logical blocks:

------------------------------------------------------------------------

## **2.1 Imports and Setup**

Required libraries:

``` python
import numpy as np
import matplotlib.pyplot as plt
import os
```

NumPy → simulation engine\
Matplotlib → plots\
OS → file management

------------------------------------------------------------------------

## **2.2 Define Utility Functions**

### **2.2.1 Periodic Laplacian**

The Laplacian drives diffusion. Periodic boundaries replicate a toroidal
microstructure model.

``` python
def laplacian_periodic(A):
    return (np.roll(A, 1, axis=0) + np.roll(A, -1, axis=0) +
            np.roll(A, 1, axis=1) + np.roll(A, -1, axis=1) - 4*A)
```

This function is used on every time step.

------------------------------------------------------------------------

## **2.3 Initialize Simulation Variables**

### **Cu concentration field**

Generated from a Gaussian distribution:

    C0 = base_C + σ * noise

Where: - `base_C` = 0.5 - `σ` = 0.10

Clipped to \[0,1\].

### **Zn mask**

Created as a random binary mask with probability = Zn fraction:

    Z = (rand < Zn_fraction).astype(int)

### **Diffusion coefficient**

Modified locally by Zn:

    D[i,j] = D0         if Z[i,j] = 0
    D[i,j] = D0 * 0.25  if Z[i,j] = 1

------------------------------------------------------------------------

## **2.4 Simulation Loop**

The physics update:

    C ← C + D * Laplacian(C) * dt

Additional steps: - clip C to \[0,1\] - compute standard deviation -
save snapshots periodically

This loop typically runs 100--300 steps.

------------------------------------------------------------------------

# 3. Output Files Generated in Step 3

Step 3 creates several important artifacts.

### **3.1 Numerical Data Files (.npz)**

For each Zn fraction, the simulation saves:

-   `C0` → Initial Cu field
-   `C_final` → Final Cu field after diffusion
-   `Z` → Zn mask
-   `stds` → Standard deviation time series
-   `times` → Matching time stamps

Example:

    step3_z0.npz
    step3_z5.npz
    step3_z10.npz

------------------------------------------------------------------------

### **3.2 Plots (.png)**

#### **3.2.1 Std vs Time**

Shows how Cu segregation decreases over annealing time.

#### **3.2.2 Cu Heatmaps**

Heatmaps at: - t = 0 - midpoint (optional) - final state

These visually confirm Zn reduces clustering.

------------------------------------------------------------------------

# 4. Example Interpretation

-   **Zn = 0%**\
    Cu clusters persist longer → slower homogenization → high std(C)

-   **Zn = 5--10%**\
    Cu mobility suppressed near Zn → faster smoothing → std(C) decreases
    faster

This matches experimental expectations from Cu-based argyrodites.

------------------------------------------------------------------------

# 5. Summary

Step 3 is where the model becomes a functional simulation tool.\
It includes:

-   full implementation of diffusion physics\
-   Zn-modified mobility\
-   metrics extraction\
-   persistent outputs\
-   essential visualization

It prepares the foundation for Step 4 (testing and experiments) and Step
5 (reporting scientific results).

------------------------------------------------------------------------

# End of Step 3 Documentation
