# Czochralski Silicon Melt — CFD Parameter Sweep Dataset

CFD data generation, case setup, mesh independence study and certification: Bertwin Kurisinkal Shine, TU Chemnitz, 2025–2026.
Produced as reference data for the Quasi project (Aditya Seshaditya, Berlin).

---

## Overview

This repository collects the certified CFD reference data from three parameter sweeps of a Czochralski (CZ) silicon melt in one place. The data is meant as ground-truth input for surrogate models and physics-informed neural networks (PINNs) developed in the Quasi project. The POD + RBF surrogate built on it is in [CZ_Surrogate_App](https://github.com/bertwinshine/CZ_Surrogate_App).

The simulation geometry and boundary conditions are based on:

> Huang, W., Cao, S., Liu, R., & Liu, D. (2025). Research on the thermal-fluid coupling in the growth process of Czochralski silicon single crystals based on an improved physics-informed neural network. *AIP Advances*, 15(10), 105202. https://doi.org/10.1063/5.0271778

All simulations were run in ANSYS Fluent 2026 R1 as 2D axisymmetric steady-state cases of the silicon melt (radius 0.30 m, height 0.15 m). Each case was exported as a CSV on the same 8181-node mesh, with the columns `r, z, u_r, u_z, u_swirl, p, T`.

---

## Parameter Sweeps

Three parameters were swept one at a time. The other two stay at the shared anchor condition.

**Shared anchor:** crystal rotation 8 rpm, crucible rotation −3 rpm, hot wall temperature 1745 K.

| Sweep | Folder | Parameter range | Cases | What the sweep shows |
|---|---|---|---:|---|
| Crystal rotation | `CZ_Crystal_Sweep/` | 4–20 rpm | 11 | the meridional flow reverses between 6 and 7 rpm |
| Hot wall temperature | `CZ_TempChange_Sweep/` | 1730–1785 K | 10 | steady flow up to 1785 K |
| Crucible rotation | `CZ_Crucible_Sweep/` | −1 to −10 rpm | 12 | velocity minimum near −4 rpm |

33 cases in total. The anchor case (crystal 8 rpm, crucible −3 rpm, 1745 K) appears in all three sweeps. The three copies are the same field; they differ only by export rounding (largest difference 5 × 10⁻⁶).

---

## Certification

Every case in this dataset passed four checks before it was included:

1. **Convergence:** all residuals settled below the prescribed threshold and showed no drift.
2. **Mesh independence:** results were confirmed stable against a refined mesh.
3. **Conservation:** mass and energy balances were verified across the domain.
4. **Boundary conditions and physics:** `certify_case.py` reads the swept value from the file name, re-derives every boundary condition from it (wall swirl, free-surface ramp, temperatures) and checks the exported field against them. It also checks that the radial and swirl velocity vanish on the axis.

Cases that did not pass all four checks were not included. To re-run check 4 on a sweep:

```
cd CZ_Crystal_Sweep
python certify_case.py *.csv
```

---

## Repository Structure

```
CZ_Crystal_Sweep/       11 cases + certify_case.py
CZ_Crucible_Sweep/      12 cases + certify_case.py
CZ_TempChange_Sweep/    10 cases + certify_case.py
```

Each sweep was first published in its own repository. Those repositories have the detailed notes: boundary-condition tables, the sweep results, certification margins and an archived release with a DOI.

| Sweep | Repository | DOI |
|---|---|---|
| Crystal rotation | [CZ_Crystal_Sweep](https://github.com/bertwinshine/CZ_Crystal_Sweep) | [10.5281/zenodo.21955299](https://doi.org/10.5281/zenodo.21955299) |
| Hot wall temperature | [CZ_Study_TempChange](https://github.com/bertwinshine/CZ_Study_TempChange) | [10.5281/zenodo.21955315](https://doi.org/10.5281/zenodo.21955315) |
| Crucible rotation | [CZ_study_Crucible_Sweep](https://github.com/bertwinshine/CZ_study_Crucible_Sweep) | [10.5281/zenodo.21955323](https://doi.org/10.5281/zenodo.21955323) |

The temperature repository also has three cases above 1785 K (1790, 1795 and 1800 K). There the melt becomes weakly unsteady, so those cases were run transient and exported as time-averaged fields. They are a different kind of data and are not part of this dataset.

---

## Intended Use

This data is meant as training and validation input for surrogate models and PINNs. All 33 cases have the same nodes in the same order, so they can be stacked directly as a snapshot matrix for POD-based surrogates, or used as reference fields for PINN residual evaluation.

If you use this dataset in your work, please cite the Huang et al. (2025) paper that the simulation setup is based on.

---

## Notes and Limitations

- Simulations were run under an ANSYS Student licence, which caps the mesh at about 512k cells. The mesh used (8181 nodes on the 2D axisymmetric domain) was confirmed mesh-independent within this limit.
- The physics model is laminar throughout the steady range. The Reynolds number of the melt convection in these cases is well below transition for this geometry.
- The free surface is imposed as a wall, not solved. No turbulence model, radiation, Marangoni convection or magnetic field.
- The dataset covers one geometry (the Huang et al. baseline). It does not cover changes in crucible geometry, melt height or material properties.
- Only one parameter moves at a time, so the data says nothing about how the three parameters interact.
- Where the anchor rotation is held fixed, it was entered as 0.838 rad/s for the crystal (8.0023 rpm) and −0.3140 rad/s for the crucible (−2.9985 rpm), 0.03 % and 0.05 % from the nominal values. All swept values are exact.
- The `certify_case.py` scripts read the swept parameter from the file name. If you rename case files, keep the naming pattern, or the checks will fail or give wrong results.

The ANSYS modelling, the sweep design and the validation are mine. The Python was written with AI assistance. I understand what it does and can walk through it.
