# Czochralski Silicon Melt — CFD Parameter Sweep Dataset

CFD data generation, case setup, mesh independence study, and certification: Bertwin Kurisinkal Shine, TU Chemnitz, 2025–2026.
Produced as reference data for the Quasi / SIMD physics-informed ML project (Aditya Seshaditya, Berlin).

---

## Overview

This repository collects the raw certified CFD reference data from three independent parameter sweeps of a Czochralski (CZ) silicon melt. The dataset was produced to serve as high-fidelity ground-truth input for surrogate models and physics-informed neural networks (PINNs) being developed under the Quasi / SIMD project.

The simulation geometry and boundary conditions are based on:

> Huang, W., Cao, S., Liu, R., & Liu, D. (2025). Research on the thermal-fluid coupling in the growth process of Czochralski silicon single crystals based on an improved physics-informed neural network. *AIP Advances*, 15(10), 105202. https://doi.org/10.1063/5.0271778

All simulations were run in ANSYS Fluent 2026 R1 as 2D axisymmetric steady-state cases of the silicon melt (radius 0.30 m, height 0.15 m). Each case was exported as a structured CSV on an 8181-node mesh with columns: `r, z, u_r, u_z, u_swirl, p, T`.

---

## Parameter Sweeps

Three axes were swept one at a time, with all other parameters held at the shared anchor condition.

**Shared anchor:** crystal rotation 8 rpm, crucible rotation −3 rpm, hot wall temperature 1745 K.

| Sweep | Parameter range | Cases | Zenodo DOI |
|---|---|---|---|
| Crystal rotation | 4–20 rpm | 11 | 10.5281/zenodo.21955299 |
| Hot wall temperature | 1730–1785 K | 10 | 10.5281/zenodo.21955315 |
| Crucible rotation | −1 to −10 rpm | 12 | 10.5281/zenodo.21955323 |

33 cases in total. The anchor case (crystal 8 rpm, crucible −3 rpm, 1745 K) appears in all three sweep repos and is bit-identical across them.

---

## Certification

Every case in this dataset passed a four-pillar certification check before inclusion:

1. **Convergence** — all residuals settled below the prescribed threshold and showed no drift.
2. **Mesh independence** — results were confirmed stable against a refined mesh.
3. **Conservation** — mass and energy balances were verified across the domain.
4. **Physical sanity** — boundary conditions (wall swirl, free-surface ramp, temperature distribution, axis symmetry) were checked against the prescribed values programmatically using `certify_case.py`.

Cases that did not pass all four checks were not included.

---

## Repository Structure

Each sweep has its own subfolder containing the certified CSV files and a dedicated README with sweep-specific notes. The `certify_case.py` script used for the certification checks is included in each sweep folder.

```
CZ_Crystal_Sweep/
CZ_Study_TempChange/
    steady/
    transient/
CZ_study_Crucible_Sweep/
```

The temperature sweep folder is split into `steady/` and `transient/` subfolders. Cases at 1785 K and below were run in steady-state mode. Cases above 1785 K showed weakly unsteady behaviour in the steady solver and were run transiently; the exported fields for those cases are time-averaged means over settled oscillation cycles.

---

## Intended Use

This data is intended as training and validation input for surrogate models and PINNs. The fields are on a consistent mesh across all 33 cases, which makes them directly usable as snapshot matrices for POD-based surrogates or as reference fields for PINN residual evaluation.

If you use this dataset in your work, please cite the Zenodo records linked above and the Huang et al. (2025) paper that the simulation setup is based on.

---

## Notes and Limitations

- Simulations were run under an ANSYS Student license, which caps the mesh at approximately 512k cells. The mesh used (8181 nodes on the 2D axisymmetric domain) was confirmed mesh-independent within this constraint.
- The physics model is laminar throughout the steady sweep range. The Reynolds number for the melt convection in these cases is well below the transition threshold for this geometry.
- The dataset covers one geometric configuration (the Huang et al. baseline geometry). It does not cover variations in crucible geometry, melt height, or material properties.
- The `certify_case.py` scripts read the swept parameter from the filename. Renaming case files without updating the filename convention will cause the certification checks to fail or produce incorrect results.
