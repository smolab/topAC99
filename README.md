# topAC99: A 99-Line MATLAB Code for Topology Optimization via the Allen--Cahn Augmented Lagrangian Method

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![SMO Educational Paper](https://img.shields.io/badge/Journal-SMO-orange.svg)](https://link.springer.com/journal/158)
[![Lab: SMOLab](https://img.shields.io/badge/Lab-SMOLab-brightgreen.svg)](https://github.com/smolab)

`topAC99.m` is a compact, self-contained, 99-line MATLAB script for structural topology optimization using the **Conserved Allen--Cahn Augmented Lagrangian Method (AC-ALM)**. 

`topAC99` lowers the entry barrier to diffuse-interface methods. It eliminates empirical density filtering (such as filter radii $r_{\min}$) and projection loops by governing spatial regularization through a physical interface width ($w$) and a double-well energy barrier ($\beta$), while capturing spontaneous topological transitions within a transparent, single-script architecture.

---

## Key Algorithmic Pillars

1. **Conforming Discrete Laplacian (Lines 21--28)**:
   Assembled in just 7 lines via orthogonal 1D Kronecker tensor sums:

   $$\mathbf{K}_{\text{lap}} = \mathbf{D}_x \otimes \mathbf{M}_y + \mathbf{M}_x \otimes \mathbf{D}_y$$

   Strictly ensures exact zero-flux Neumann boundary conservation ($\sum_j K_{\text{lap}, ij} = 0$), eliminating numerical boundary phase drift without auxiliary projection loops.
2. **Matrix-Free Convolution Architecture (Lines 35, 43)**:
   Employs native 2D spatial convolution (`conv2`) with a $2 \times 2$ kernel ($1/4$) for physical element density mapping and adjoint sensitivity backpropagation, completely bypassing large sparse interpolation matrices.
3. **Closed-Form 1D-Dual MMA Solver (Lines 58--85)**:
   Operates on a guaranteed finite dual bracket $\lambda \in [0, 10^6]$ in 28 lines. The solver converges via bisection in 25 iterations to machine precision, naturally handling signed phase-separation sensitivities while strictly preserving physical volume feasibility.

---

## Quick Start

### Function Signature
```matlab
[c, x, loop, h_norm] = topAC99(nelx, nely, volfrac, penal, w, beta)
```

### Standard Half-MBB Benchmark
Run the standard half-MBB beam ($120 \times 80$ mesh, volume fraction $V_0 = 0.40$):
```matlab
topAC99(120, 80, 0.40, 3.0, 2.1, 0.065);
```

**Expected Console Output**:
```text
It.: 81  Comp.: 59.163  ||h||_*: 5.62e-03  ch.: 0.0099
```
* Optimization terminates cleanly in **81 iterations** at Sigmund's classic threshold $\text{ch} < 0.01$ without numerical lag.
* Optimized compliance: $c = 59.163$.
* Exact volume error: $|V - V_0| = 5.1 \times 10^{-6}$ ($< 0.0006\%$ error).
* Phase equilibrium dual residual: $\|\mathbf{h}\|_* = 5.62 \times 10^{-3}$.

For target volume fraction $V_0 = 0.30$:
```matlab
topAC99(120, 80, 0.30, 3.0, 2.1, 0.065);
```
**Expected Output**: Converges in **80 iterations** ($c = 77.435$, $\|\mathbf{h}\|_* = 6.90 \times 10^{-3}$).

---

## Parameter Guidelines

| Parameter | MATLAB Variable | Baseline | Recommended Range | Description & Physical Role |
| :--- | :--- | :--- | :--- | :--- |
| **Mesh Dimensions** | `nelx`, `nely` | `120, 80` | Any grid size | Cartesian element count along horizontal and vertical axes. |
| **Volume Fraction** | `volfrac` | `0.40` | $0.10 \sim 0.60$ | Prescribed target material volume fraction $V_0$. |
| **Penalization Power** | `penal` | `3.0` | $3.0 \sim 4.0$ | SIMP stiffness interpolation power ($p$). |
| **Interface Width** | `w` | `2.1` | $1.5 \sim 4.2$ | Observable diffuse transition width $w = \xi_{0.1\to 0.9} = \sqrt{2}\ln(9)\ell_{\text{AC}}$ (in elements). |
| **Energy Barrier** | `beta` | `0.065` | $0.03 \sim 0.12$ | Double-well potential scale $W(\phi) = \beta\phi^2(1-\phi)^2$; sets phase-segregation driving force. |

### Output Arguments
* `c`: Optimized structural mean compliance ($\mathbf{F}^\top \mathbf{u}$).
* `x`: Optimized physical element density matrix ($n_y \times n_x$).
* `loop`: Total iteration count to convergence.
* `h_norm`: Mass-weighted continuous dual RMS residual norm $\|\mathbf{h}\|_* = \sqrt{\frac{1}{|\Omega|} \sum_i h_i^2 / m_i}$.

---

## Modular Extensions

`topAC99.m` easily adapts to diverse boundary conditions and passive design regions:

### 1. Clamped Cantilever Beam with Center Tip Load
Modify Lines 18--19 to apply a fixed left wall and centered tip point load:
```matlab
% Cantilever Boundary Condition & Load (replace Lines 18--19)
F = sparse(2*(nelx*(nely+1) + round(nely/2)+1), 1, -1, 2*(nely+1)*(nelx+1), 1);
fixeddofs = 1:2*(nely+1);
```
Invoking `topAC99(120, 80, 0.40, 3.0, 2.1, 0.065)` converges in **77 iterations** ($c = 43.992$).

### 2. Highway Deck Arch Bridge with Passive Solid Deck
Model a $120 \times 40$ bridge with two pinned bedrock abutments, distributed top roadway deck load, and a 2-element non-optimizable passive solid roadway slab:

First, insert the following lines immediately after Line 72 to define the passive solid deck region:
```matlab
% Passive Solid Roadway Deck (insert immediately after Line 72)
[node_X, node_Y] = meshgrid(0:nelx, 0:nely);
a_lim(node_Y <= 2) = 1.0; b_lim(node_Y <= 2) = 1.0;
```
Second, replace Lines 18--19 with the following to apply the bridge boundary conditions and distributed deck load:
```matlab
% Highway Deck Arch Bridge Load & supports (replace Lines 18--19)
deck_nodes = 1:( nely +1) :( nelx +1) *( nely +1);
f_deck = [0.5; ones ( nelx -1 ,1) ; 0.5] / nelx;
F = sparse (2* deck_nodes , 1 , - f_deck , 2*( nely +1) *( nelx +1) , 1);
fixeddofs = [2*( nely +1) +( -1:0) , 2*( nelx +1) *( nely +1) +( -1:0)];
```
Invoking `topAC99(120, 40, 0.40, 3.0, 2.1, 0.065)` converges in **88 iterations** ($c = 6.447$).

---

## Citation

If you use `topAC99` in your research, teaching, or academic projects, please cite:

```bibtex
@article{Du2026topAC99,
  author  = {Run Du},
  title   = {{topAC99}: A 99-Line {MATLAB} Code for Topology Optimization via the {Allen--Cahn} Augmented {Lagrangian} Method},
  journal = {Structural and Multidisciplinary Optimization},
  year    = {2026},
  note    = {Educational Paper, under review}
}
```

---

## About the Author & Laboratory

* **Author**: Run Du ([rdu@swjtu.edu.cn](mailto:rdu@swjtu.edu.cn))
* **Laboratory**: **SMOLab** (Structural & Multidisciplinary Optimization Laboratory)
* **Affiliations**: 
  * School of Mechanical Engineering, Southwest Jiaotong University, Chengdu, China
  * State Key Laboratory of Bridge Intelligent and Green Construction, Southwest Jiaotong University
* **Official Repository**: [https://github.com/smolab/topAC99](https://github.com/smolab/topAC99)

---

## License

This software is released under the **GNU Affero General Public License v3.0 (AGPLv3)**. See the `topAC99.m` bottom, `LICENSE` or the [GNU website](https://www.gnu.org/licenses/) for details.
