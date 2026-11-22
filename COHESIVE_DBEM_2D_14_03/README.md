# Multizone Cohesive DBEM 2D (`COHESIVE_DBEM_2D`)

Modern Fortran computational code (**OpenMP Multithread + LAPACK/BLAS**) with a Python preprocessing and postprocessing pipeline (**Gmsh + `meshio` + `matplotlib`**) for linear elastic analysis and interfacial fracture/debonding propagation in adhesively bonded joints and 2D composite structures, using the **Multizone Cohesive Dual Boundary Element Method (Multizone Cohesive DBEM 2D)**.

This repository and its formulation are directly associated with the scientific paper:

> **Comparative analysis of adhesive debonding in three-point bending tests using boundary element methods: Coupled criterion versus interface damage model**  
> *Sérgio Gustavo Ferreira Cordeiro, Thiago V. Birro, René Quispe Rodriguez, Eder Lima de Albuquerque, Lucas Silveira Campos, Jairo Francisco Useche Vivero, Éric Paroissien, Maëlenn Aufray, Frédéric Lachaud.*  
> **Engineering Analysis with Boundary Elements (EABE)**, vol. 180, 106474, 2025.  
> Reference file: [ArtigoEABE2025.pdf](../ArtigoEABE2025.pdf)

---

## 1. Theoretical Formulation Overview

### 1.1 Multizone Dual Boundary Element Method (Multizone DBEM)
Unlike the Finite Element Method (FEM), which requires volumetric discretization of the entire domain and suffers from severe mesh distortion in thin adhesive layers, the BEM discretizes only the boundaries and interfaces of the bodies.

For heterogeneous systems (e.g., a metallic substrate bonded to a polymeric adhesive block), the domain is decomposed into homogeneous, isotropic subregions $\Omega_1, \Omega_2, \dots, \Omega_N$. Within each subregion, the Displacement Boundary Integral Equation (DBIE) and the Traction Boundary Integral Equation (TBIE) are written as:

$$c_{ij}(\mathbf{x}') u_j(\mathbf{x}') + \text{CPV} \int_{\Gamma} T_{ij}^*(\mathbf{x}', \mathbf{x}) u_j(\mathbf{x}) \, \mathrm{d}\Gamma = \int_{\Gamma} U_{ij}^*(\mathbf{x}', \mathbf{x}) t_j(\mathbf{x}) \, \mathrm{d}\Gamma$$

$$c_{ij}(\mathbf{x}') t_j(\mathbf{x}') + n_k(\mathbf{x}') \, \text{HFP} \int_{\Gamma} S_{kij}^*(\mathbf{x}', \mathbf{x}) u_j(\mathbf{x}) \, \mathrm{d}\Gamma = n_k(\mathbf{x}') \, \text{CPV} \int_{\Gamma} D_{kij}^*(\mathbf{x}', \mathbf{x}) t_j(\mathbf{x}) \, \mathrm{d}\Gamma$$

where:
- $\mathbf{x}'$ and $\mathbf{x}$ are the source (collocation) and field points;
- $U_{ij}^*, T_{ij}^*, D_{kij}^*, S_{kij}^*$ are the Kelvin fundamental tensors for 2D elasticity (plane stress `PT` or plane strain `PD`);
- $\text{CPV}$ and $\text{HFP}$ denote integrals evaluated in the Cauchy Principal Value and Hadamard Finite Part senses, respectively.

### 1.2 Discontinuous and Semi-Discontinuous Elements
To avoid ambiguities in the unit outward normal vector at geometric corners and boundary/interface transitions:
- Double nodes are placed at corners and interface endpoints;
- On discontinuous elements (`D`, `LD`, `RD`), the **collocation points** are located inside the parametric element ($\xi = \pm 1/3$ for linear elements), ensuring smooth and accurate evaluation of interfacial tractions.

### 1.3 Interfacial Coupling and Progressive Interface Damage Model (PID)
Interfaces between subregions are divided into paired surfaces:
- Surface **$I^-$ (`Im`)**: associated with subregion 1 (e.g., substrate);
- Surface **$I^+$ (`Ip`)**: associated with subregion 2 (e.g., adhesive).

The solver supports three analysis modes configured in `4_ANALYSIS_PARAMETERS.txt`:
1. **Linear Elastic Analysis (`LE`)**: evaluates normal peel ($\sigma$) and shear ($\tau$) stresses along the undamaged weak interface with penalized stiffness ($k_{n0}, k_{t0}$), used in the **Coupled Criterion (CC)** inverse analysis.
2. **Symmetric Nonlinear Cohesive Propagation (`CP`)**: employs a mixed-mode bilinear traction-separation law (Mode I: $K_{n0}, \Phi_n, \sigma_{\max}$; Mode II: $K_{t0}, \Phi_t, \tau_{\max}$) with a scalar damage variable $D \in [0, 1]$ solved incrementally via **Newton-Raphson** continuation with displacement or relative-opening control.
3. **Asymmetric Nonlinear Cohesive Propagation (`CA`)**: reproduces the asymmetric debonding observed experimentally in the 3PBT (Section 6.3.2 of the EABE 2025 paper) by assigning $\sigma_c = 10^{50}\text{ MPa}$ to the right corner of the interface so that a single crack nucleates and propagates exclusively from the left corner.

---

## 2. Modular Architecture and Performance Optimizations (OpenMP)

The Fortran solver has been reviewed and optimized for high performance while preserving the mathematical clarity of the formulation:

| Module | File | Responsibility and Optimizations |
|---|---|---|
| **Main Program** | `MAIN_COHESIVE_DBEM_2D.f90` | OpenMP environment initialization and wall-clock timing via `omp_get_wtime()`. |
| **Global Data** | `MODULE_DATA.f90` | Double-precision (`real64`) global structures, $G$ and $H$ matrices, connectivities, and material properties. |
| **Input Data** | `MODULE_INPUT_DATA.f90` | Reading of the 5 `.txt` input files, Gauss-Legendre quadrature, and leak-free topological mapping. |
| **Data Treatment** | `MODULE_DATA_TREATMENT.f90` | Discontinuous collocation point computation and parallel (`!$OMP PARALLEL DO`) pairing of $I^- \leftrightarrow I^+$ cohesive nodes. |
| **Geometric Mappings** | `MODULE_GEOMETRIC_MAPPINGS.f90` | Lagrange geometric interpolation, Jacobian (`DSQRT`), and outward normal vector. |
| **Shape Functions** | `MODULE_APPROXIMATION_FUNCTIONS.f90` | Basis functions and derivatives at nodal and collocation points. |
| **Numerical Integration** | `MODULE_NUMERICAL_INTEGRATION.f90` | Regular and singular (CPV/HFP) integration for DBIE and TBIE using stack-allocated automatic arrays (eliminating $\sim 400{,}000$ dynamic heap allocations) and hoisted elastic constants outside the Gauss loop. |
| **Analysis & Assembly** | `MODULE_ANALYSIS.f90` | Assembly of $H$ and $G$ matrices parallelized with **`!$OMP PARALLEL DO SCHEDULE(DYNAMIC)`** and an $O(1)$ node-to-element lookup table (`NODE_ELEM`). |
| **Linear / Nonlinear Solver** | `MODULE_SOLVER.f90` | Parallel assembly of the coupled system matrix, sparse evaluation of prescribed boundary conditions, linear solution via LAPACK (`DGESV`), smart LU factorization caching (`DGETRF`/`DGETRS`) for nonlinear Newton-Raphson iterations, and step-by-step displacement export (`12_OUTPUT_NODAL_DISPLACEMENTS_STEPS.txt`). |
| **Postprocessing** | `MODULE_POSTPROCESSING.f90` | OpenMP-parallelized evaluation of displacements and stresses at internal and boundary points, exporting `11_OUTPUT_BOUNDARY_STRUCTURAL_RESPONSE.txt` and portable `OUTPUT.vtk`. |

### Achieved Speedup
On the EABE 2025 benchmark ($N_{\text{nodes}} = 445$ boundary nodes, $N_{\text{elem}} = 445$ elements, $N_{\text{Gauss}} = 200$ quadrature points per element, 8-thread CPU):
- **Original version (`-O0 -fbounds-check`, sequential)**: $\approx 38.05\text{ s}$ per simulation.
- **Compiler optimization (`-O3`, sequential)**: $\approx 11.75\text{ s}$ per simulation ($3.2\times$ faster).
- **OpenMP parallelization + stack arrays + `-O3`**: **$\approx 2.97\text{ s}$ per simulation** (**$12.8\times$ faster**, preserving identical numerical results to the 9th decimal place).

---

## 3. Interactive Jupyter Notebooks (Local & Google Colab)

The project includes three complete Jupyter Notebooks covering both linear elastic and nonlinear incremental progressive interface damage (PID) analyses:

1. **`Run_cohesive_dbem_2d.ipynb`** (*Linear Elastic Analysis — Local / Cross-Platform*):
   - Runs the entire workflow from parametric `.geo` creation in Gmsh, `.msh` mesh generation, `meshio` preprocessing, Fortran compilation, and linear elastic solving for $h_s = 1.08\text{ mm}$ and $h_s = 1.62\text{ mm}$, reproducing **Fig. 10** and **Figs. 12, 13, 14, and 15** of the paper.
2. **`Run_cohesive_dbem_2d_on_Colab.ipynb`** (*Linear Elastic Analysis — Google Colab Linux*):
   - Specifically designed for execution on **Google Colab (Ubuntu Linux)** using direct shell commands (`!apt-get`, `!make`, `!gmsh`, `!./MAIN_COHESIVE_DBEM_2D`).
3. **`Run_cohesive_dbem_2d_nlinear.ipynb`** (*Nonlinear Progressive Interface Damage Analysis — Step by Step*):
   - Details step by step the theory and implementation of displacement-controlled incremental nonlinear analysis (Newton-Raphson with a mixed-mode bilinear cohesive law) for asymmetric debonding (`ANALYSIS_TYPE = 'CA'`), reproducing **Fig. 26** (*Force–displacement results for asymmetric debonding: $h_s = 1.08\text{ mm}$*), crack length evolution (analogous to **Fig. 30**), and the synchronized `.mp4` animation.

---

## 4. EABE 2025 Benchmark (Three-Point Bending Test — Fig. 10 & Fig. 26)

### 4.1 Geometric and Physical Parameters
- **Substrate ($\Omega_1$)**: length $L_s = 50.0\text{ mm}$, thickness $h_s = 1.08\text{ mm}$ or $h_s = 1.62\text{ mm}$, $E_1 = 136{,}000\text{ MPa}$ ($2 \times 68{,}000\text{ MPa}$ scaled by the width ratio $b_1/b_2$), $\nu_1 = 0.33$.
- **Adhesive ($\Omega_2$)**: bonded overlap length $L_a = 25.0\text{ mm}$ ($x \in [-12.5, +12.5]\text{ mm}$), thickness $h_a = 4.0\text{ mm}$, $E_2 = 3{,}000\text{ MPa}$, $\nu_2 = 0.35$.
- **Supports and Loading**: span between supports $d_{\text{supp}} = 35.90\text{ mm}$ (left rolling support $u_y = 0$, right fixed support $u_x = u_y = 0$, width $0.12\text{ mm}$) and central top loading zone ($w_{\text{load}} = 2.0\text{ mm}$):
  - For $h_s = 1.08\text{ mm}$: $U_c = 0.576\text{ mm}$ ($F_c = 178\text{ N}$, $\sigma_c = 44.654\text{ MPa}$, $\mathcal{G}_c = 0.520\text{ N/mm}$, $k_{n0} = k_{t0} = 2000\text{ MPa/mm}$);
  - For $h_s = 1.62\text{ mm}$: $U_c = 0.6093\text{ mm}$ ($F_c = 326\text{ N}$, $\sigma_c = 45.906\text{ MPa}$, $\mathcal{G}_c = 0.5371\text{ N/mm}$, $k_{n0} = k_{t0} = 948\text{ MPa/mm}$).

### 4.2 Automatically Generated Figures & Animation
- `mesh_and_boundary_conditions.png`: Reproduction of **Fig. 10** (mesh and boundary conditions for linear analysis, $h_s = 1.08\text{ mm}$).
- `mesh_and_boundary_conditions_nlinear.png`: Mesh, boundary conditions, and control node for displacement-controlled nonlinear analysis ($h_s = 1.08\text{ mm}$).
- `fig12_peel_hs108.png`: **Fig. 12** — *Peel stress over the undamaged interface: $h_s = 1.08\text{ mm}$* (peak $\approx 43.45\text{ MPa}$ at the edges $x = \pm 12.5\text{ mm}$).
- `fig13_shear_hs108.png`: **Fig. 13** — *Shear stress over the undamaged interface: $h_s = 1.08\text{ mm}$* (antisymmetric extrema $\approx \mp 13.96\text{ MPa}$).
- `fig14_peel_hs162.png`: **Fig. 14** — *Comparison of peel Stress over the undamaged interface: $h_s = 1.62\text{ mm}$* (peak $\approx 46.15\text{ MPa}$).
- `fig15_shear_hs162.png`: **Fig. 15** — *Comparison of shear stress over the undamaged interface: $h_s = 1.62\text{ mm}$* (extrema $\approx \mp 15.09\text{ MPa}$).
- `fig26_force_displacement_hs108.png`: **Fig. 26** — *Force–displacement results for asymmetric debonding: $h_s = 1.08\text{ mm}$* (nonlinear **2D BEM - PID** curve).
- `fig26_cohesive_damage_evolution_hs108.png`: Asymmetric crack length evolution $a(u)$ (analogous to **Fig. 30**) and normalized cohesive stiffness profile $K_{nn}/k_{n0}$.
- `debonding_animation_hs108.mp4` (generated by `animate_debonding.py`): Synchronized `.mp4` animation showing the 2D DBEM mesh deforming on a $1:1$ scale (**Fig. 28**) alongside the *Force × Central displacement* graph (**Fig. 26**).

---

## 5. Command-Line Build and Execution

### Windows (MSYS2 MinGW-w64 UCRT64)
```powershell
$env:PATH = "C:\msys64\ucrt64\bin;C:\msys64\usr\bin;" + $env:PATH
pip install -r requirements.txt
mingw32-make clean
mingw32-make
.\MAIN_COHESIVE_DBEM_2D.exe
```

### Linux / Google Colab
```bash
sudo apt-get update && sudo apt-get install -y gfortran liblapack-dev libblas-dev gmsh
pip install -r requirements.txt
make clean && make
./MAIN_COHESIVE_DBEM_2D
```
