
# Development of Stable Cu-Based Argyrodites for Highly Efficient Energy Conversions

**Researcher:** Sagarika Sharma  
**Project No.:** 10 
**Programming Method Used:** Python  

---



## 1. Background and Key Concepts

### 1.1 Cu-based Argyrodites
- Argyrodites are crystalline materials with general formula **Cu₆PS₅X** (X = Cl, Br, I).  
- They exhibit **superionic conduction**, meaning Cu⁺ ions are highly mobile within the lattice.  
- Applications include **thermoelectrics** and **solid-state ionic conductors**.

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

## 2. Computational Model (Simple 2D Toy Model)

**Goal:** Demonstrate the effect of Zn on Cu segregation visually and quantitatively.

### 2.1 Model Setup
- 2D grid (e.g., 50×50 cells) represents material microstructure  
- Each cell contains:
  - \(C_{i,j}\): local Cu fraction (0–1)  
  - \(Z_{i,j}\): Zn presence (0 or 1)  

### 2.2 Update Rule
- **Diffusion-like smoothing** of Cu:
  \[
  C^{t+1}_{i,j} = C^t_{i,j} + D_{i,j} \cdot (\text{Laplacian of } C)
  \]
- Local diffusion coefficient:
  \[
  D_{i,j} =
    \begin{cases}
      D_0 & \text{if no Zn} \\
      D_0 \cdot 0.25 & \text{if Zn present}
    \end{cases}
  \]

---

## 3. Work Plan with Milestones and Fine-Grained Tasks

| Milestone | Task | Description | Estimated Time |
|-----------|------|-------------|----------------|
| **0: Setup** | 0.1 | Install Python and libraries (`numpy`, `matplotlib`) | 20–30 min |
|  | 0.2 | Create project folder and save script | 5–10 min |
| **1: Literature Review & Parameter Selection** | 1.1 | Read key papers on Cu segregation and Zn substitution | 4–6 hours |
|  | 1.2 | Decide grid size, diffusion coefficient, Zn fraction | 1–2 hours |
|  | 1.3 | Document assumptions and parameter table | 1 hour |
| **2: Model Design** | 2.1 | Write pseudocode for diffusion update and metrics | 2–3 hours |
|  | 2.2 | Decide outputs: plots, CSV files, snapshots | 1 hour |
| **3: Python Prototype Implementation** | 3.1 | Implement diffusion update and Zn mask | 3–5 hours |
|  | 3.2 | Add metrics calculation (std, optional CV) | 2–3 hours |
|  | 3.3 | Implement simulation loop and collect results | 2–3 hours |
|  | 3.4 | Generate plots (time series and maps) | 2 hours |
| **4: Testing and Experiments** | 4.1 | Run simulation with varying Zn fractions (0, 0.05, 0.1) | 30–60 min |
|  | 4.2 | Inspect and validate plots | 20–30 min |
| **5: Reporting Results** | 5.1 | Screenshot figures and add captions | 20–30 min |
|  | 5.2 | Write short report paragraph explaining results | 20–30 min |
| **6: Optional Extensions** | 6.1 | Multi-run statistics and confidence intervals | 2–3 hours |
|  | 6.2 | Cluster size or spatial autocorrelation analysis | 2–3 hours |

**Total estimated time:** ~2–3 hours for basic prototype; ~1–2 days if including optional statistics/plots

---

## 4. Schedule Proposal

**Day 1 (2–3 hours)**  
- Install Python, save script  
- Run prototype simulation  
- Generate std(Cu) plot and Cu maps  

**Day 2 (1–2 hours)**  
- Test different Zn fractions  
- Screenshot figures  
- Write short report paragraph  

**Optional Day 3**  
- Run multiple simulations for statistics  
- Compute confidence intervals and cluster analysis  
- Polish figures for report or presentation  

---

## 5. Expected Results

- **Time series plot:** std(Cu) decreases faster with Zn → reduced segregation  
- **Cu maps:** more uniform Cu distribution in Zn-substituted grid  
- Supports the hypothesis that **Zn stabilizes Cu distribution during annealing**  

---

## 6. Sample Report Paragraph

> We simulated Cu redistribution on a 2D lattice during annealing, comparing two cases: no Zn and 5% Zn substitution. The standard deviation of Cu concentration decreased faster in the Zn-substituted lattice, indicating lower segregation. Visual inspection of concentration maps confirmed that Zn reduced clustering, producing a more uniform distribution. This simple model demonstrates that Zn reduces Cu mobility and stabilizes the microstructure, explaining improved reproducibility of transport properties in Zn-doped argyrodites.
