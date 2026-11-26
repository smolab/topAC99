# `topAC99`: A compact MATLAB implementation of filter-free Allen--Cahn topology optimization with matrix-free field transfer

[![MATLAB](https://img.shields.io/badge/MATLAB-R2018b+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![GNU Octave](https://img.shields.io/badge/GNU%20Octave-7.0+-0790C0.svg)](https://octave.org/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)
[![Lab: SMOLab](https://img.shields.io/badge/Lab-SMOLab-brightgreen.svg)](https://github.com/smolab)

`topAC99.m` is a self-contained, toolbox-free 99-line MATLAB kernel for structural topology optimization using **Allen--Cahn Residual-Constrained Topology Optimization (ACTO)** via a mass-weighted **Augmented Lagrangian Method (AC-ALM)**.

Operating out-of-the-box in both base MATLAB (R2018b+) and GNU Octave (7.0+) using only built-in array and sparse primitives (`conv2`, `sparse`, `spdiags`, `kron`, `mldivide`), `topAC99` eliminates empirical density filtering (such as filter radii $r_{\min}$), auxiliary Helmholtz PDE solves, and external optimization packages. Spatial regularity is governed intrinsically by a physical diffuse-interface width $w$ and a double-well barrier $\beta$.

---

## Key Algorithmic Components

1. **Lumped Tensor-Product Laplacian & Mass Operators (Lines 21--28)**:
   Assembled in 7 lines via orthogonal 1D Kronecker tensor sums with trapezoidal nodal quadrature weights $\mathbf{M}$:

   $$\mathbf{K}_{\mathrm{lap}} = \mathbf{D}_x \otimes \mathbf{M}_y + \mathbf{M}_x \otimes \mathbf{D}_y$$

   This stencil enforces exact zero-flux Neumann boundary conditions ($\mathbf{1}^\top \mathbf{K}_{\mathrm{lap}} = \mathbf{0}^\top$) with width-invariant $\mathcal{O}(N_{\mathrm{n}})$ sparse storage ($6.28\text{ MB}$ on a $300 \times 200$ grid regardless of interface width $w$).

2. **Matrix-Free Dual-Grid Field Transfer via `conv2` (Lines 35, 43) & Staggered AC-ALM (Lines 45--56, 93--96)**:
   Primal four-node centroid averaging (`conv2(phi, [1 1; 1 1]/4, 'valid')`) and exact adjoint sensitivity backpropagation (`conv2(dc, [1 1; 1 1]/4, 'full')`) bypass explicit sparse node--element transfer matrices, paired with a staggered continuation-weighted ALM update reusing the evaluated residual $\mathbf{h}$.

3. **Inlined 18-Line 1D-Dual Bisection MMA Solver (Lines 67--85; Asymptotes Lines 58--65)**:
   A closed-form separable convex subproblem solver parameterized by 25 bisection steps on $\lambda \in [0, 10^6]$ with a conservative volume surrogate and an analytical lower-bound safeguard (`phi_feas = a_lim`), natively accommodating signed sensitivities without external optimizer packages.

---

## Quick Start

### Function Signature
```matlab
[c, x, loop, h_norm] = topAC99(nelx, nely, volfrac, penal, w, beta)
```

### Standard Half-MBB Benchmarks
Run the $120 \times 80$ half-MBB beam at $V_0 = 0.30$ and $V_0 = 0.40$:
```matlab
[c, x, loop, h_norm] = topAC99(120, 80, 0.30, 3.0, 2.1, 0.065);
[c, x, loop, h_norm] = topAC99(120, 80, 0.40, 3.0, 2.1, 0.065);
```

**Expected Reference Outputs** (single-thread MATLAB R2023a, dual stopping gate $\text{loop} \ge 30$, $\text{ch} < 0.01$, $\|\mathbf{h}\|_* \le 0.28\beta$):
* $V_0 = 0.30$: converges in **80 iterations** with $c = 77.435$, thresholded binary compliance $c_{0/1} = 70.616$, non-discreteness $M_{\mathrm{nd}} = 11.16\%$, dual residual $\|\mathbf{h}\|_* = 6.92 \times 10^{-3}$, and volume error $|V - V_0| = 3.09 \times 10^{-6}$.
* $V_0 = 0.40$: converges in **81 iterations** with $c = 59.163$, $c_{0/1} = 55.416$, $M_{\mathrm{nd}} = 11.64\%$, $\|\mathbf{h}\|_* = 5.62 \times 10^{-3}$, and $|V - V_0| = 5.10 \times 10^{-6}$.
* Multi-resolution benchmarks up to $300 \times 200$ ($121{,}002$ elastic DOFs) confirm that Allen--Cahn terms add only $0.7\%$--$1.2\%$ per-step overhead onto FEA-dominated iterations while converging in $78$--$92$ iterations on fine grids ($120 \times 80$ to $300 \times 200$).

---

## Parameter Guidelines

| Parameter | MATLAB Variable | Baseline | Recommended Range | Description & Physical Role |
| :--- | :--- | :--- | :--- | :--- |
| **Mesh Dimensions** | `nelx`, `nely` | `120, 80` | Any structured Q4 grid | Horizontal and vertical element counts ($n_x \times n_y$). |
| **Volume Fraction** | `volfrac` | `0.30` / `0.40` | $0.10 \sim 0.60$ | Target material volume fraction $V_0$. |
| **Penalization Power** | `penal` | `3.0` | $3.0 \sim 4.0$ | SIMP stiffness interpolation exponent $p$. |
| **Interface Width** | `w` | `2.1` | $2.0 \sim 4.2$ | Diffuse transition width $w = \sqrt{2}\ln(9)\ell_{\mathrm{AC}}$ (in elements; $w \ge 2.0$ resolves the $\tanh$ profile across multiple elements). |
| **Energy Barrier** | `beta` | `0.065` | Problem-scaled | Double-well potential coefficient $W(\phi) = \beta\phi^2(1-\phi)^2$. |

* **Order-of-Magnitude $\beta$ Calibration**: When changing objective scales or boundary conditions, $\beta$ should roughly balance the initial penalty restoring gradient against the mechanical sensitivity norm:
  $$\beta \sim \mathcal{O}\!\left(\sqrt{\frac{\|\nabla_{\boldsymbol{\phi}} c\|_{\mathrm{rms}}}{w_h \rho_{h,0}}}\right)$$
  For example, switching from half-MBB compliance ($\beta = 0.065$) to the unit-force compliant inverter ($\sim 10^3\times$ smaller sensitivity norm) scales $\beta$ by $\sim 1/\sqrt{10^3}$ to $\beta = 0.002$.
* **Ultra-Low Volume Fractions ($V_0 \le 0.15$)**: Reducing the initial trust radius `move` on Line 32 to $0.2V_0$--$0.3V_0$ prevents early step oscillations from severing nascent load paths.

### Output Arguments
* `c`: Structural objective value at convergence (mean compliance $\mathbf{F}^\top \mathbf{u}$).
* `x`: Physical element density matrix ($n_y \times n_x$).
* `loop`: Total iteration count to convergence.
* `h_norm`: Mass-weighted continuous dual RMS residual norm $\|\mathbf{h}\|_* = \sqrt{\frac{1}{|\Omega|} \sum_i h_i^2 / m_i}$.

---

## Minimal-Line Extensions

`topAC99.m` adapts directly to non-self-adjoint compliant mechanisms and passive structural subdomains without modifying the Kronecker Laplacian, `conv2` field transfers, AC-ALM loop, or 18-line 1D-dual MMA solver:

### 1. Compliant Mechanism: Symmetric Force Inverter
Apply the edits below from bottom to top so the earlier line numbers in `topAC99.m` stay valid (7 kernel lines replaced by 9 lines):

1. First, replace **Lines 38--41** with the port springs, two-column solve, and mutual-energy evaluation:
```matlab
% Springs, Two-Column Solve & Mutual Energy (Replaces Lines 38--41)
K = sparse([iK; din; dout], [jK; din; dout], [sK; 0.1; 0.1]); K = (K + K')/2;
U(freedofs,:) = K(freedofs, freedofs) \ F(freedofs,:);
U1 = U(:,1); U2 = U(:,2);
ce = reshape(sum((U1(edofMat)*KE).*U2(edofMat), 2), nely, nelx);
c = U(dout, 1);
```
2. Next, replace **Lines 17--19** with the port definitions, loads, and symmetry boundary conditions:
```matlab
% Force Inverter Ports, Loads & Symmetry BCs (Replaces Lines 17--19)
din = 1; dout = 2*nelx*(nely+1) + 1;
U = zeros(2*(nely+1)*(nelx+1), 2);
F = sparse([din, dout], [1, 2], [1, 1], 2*(nely+1)*(nelx+1), 2);
fixeddofs = [2*(nely+1)+(-3:0), 2:2*(nely+1):2*(nelx+1)*(nely+1)];
```
Running `topAC99(120, 60, 0.30, 3.0, 2.1, 0.002)` on the lower symmetric half-grid converges in **93 iterations** to $u_{\mathrm{out}} = -1.0180$ ($u_{\mathrm{out},0/1} = -1.0606$, geometric advantage $\mathrm{GA} = 0.115$, $\|\mathbf{h}\|_* = 3.69 \times 10^{-4}$, $|V - V_0| = 1.17 \times 10^{-4}$). Homogeneous zero-flux symmetry across row 1 holds automatically via `Klap`.

### 2. Highway Deck Arch Bridge with Passive Solid Deck
Model a $120 \times 40$ highway deck arch bridge ($V_0 = 0.40$) with pinned bedrock abutments, a uniform downward deck load, and a continuous two-element passive solid roadway deck (nodal rows `1:3`). Apply the edits from bottom to top to preserve line numbers:

1. First, insert **after Line 72** to clamp the passive deck bounds:
```matlab
a_lim(1:3, :) = 1.0; b_lim(1:3, :) = 1.0;      % insert after Line 72
```
2. Second, insert **after Line 30** to initialize the passive deck and history buffers:
```matlab
phi(1:3, :) = 1.0; p_old1 = phi; p_old2 = phi; % insert after Line 30
```
3. Third, replace **Lines 18--19** with the distributed deck load and pinned bedrock supports:
```matlab
% Distributed Deck Load & Bedrock Pins (Replaces Lines 18--19)
deck_nodes = 1:(nely+1):(nelx+1)*(nely+1);
f_deck = [0.5; ones(nelx-1,1); 0.5] / nelx;
F = sparse(2*deck_nodes, 1, -f_deck, 2*(nely+1)*(nelx+1), 1);
fixeddofs = [2*(nely+1)+(-1:0), 2*(nelx+1)*(nely+1)+(-1:0)];
```
Running `topAC99(120, 40, 0.40, 3.0, 2.1, 0.065)` converges in **87 iterations** to $c = 6.501$ ($c_{0/1} = 5.938$, $|V - V_0| = 6.36 \times 10^{-5}$).

---

## Citation

If you use `topAC99` in your research or teaching, please cite both this software implementation and the companion continuum ACTO formulation:

```bibtex
@article{DuRun2026topAC99,
  author  = {Run Du},
  title   = {{topAC99}: A compact {MATLAB} implementation of filter-free {A}llen--{C}ahn topology optimization with matrix-free field transfer},
  journal = {arXiv preprint arXiv:2611.31885},
  year    = {2026},
  doi     = {10.48550/arXiv.2611.31885}
  note    = {Submitted manuscript; code available at \url{https://github.com/smolab/topAC99}}
}

@article{DuRun2026topac,
  author  = {Run Du and Xueke Xiao and Wenming Cheng and Zeqiang Zhang},
  title   = {{A}llen--{C}ahn residual-constrained topology optimization via an augmented {L}agrangian method},
  journal = {arXiv preprint arXiv:2611.02375},
  year    = {2026},
  doi     = {10.48550/arXiv.2611.02375}
}
```

---

## Author & Laboratory

* **Author**: Run Du ([rdu@swjtu.edu.cn](mailto:rdu@swjtu.edu.cn))
* **Laboratory**: **SMOLab** (Structural & Multidisciplinary Optimization Laboratory)
* **Affiliations**:
  * School of Mechanical Engineering, Southwest Jiaotong University, Chengdu, China
  * Technology and Equipment of Rail Transit Operation and Maintenance Key Laboratory of Sichuan Province
  * State Key Laboratory of Bridge Intelligent and Green Construction
* **Repository**: [https://github.com/smolab/topAC99](https://github.com/smolab/topAC99)

---

## License

Released under the **GNU Affero General Public License v3.0 (AGPLv3)**. See [LICENSE](LICENSE) or the footer of [`topAC99.m`](topAC99.m) for details.
