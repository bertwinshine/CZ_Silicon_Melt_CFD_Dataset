# Czochralski silicon melt: CFD parameter sweep dataset

> **Superseded by [CZ-Dataset-Doe](https://github.com/bertwinshine/CZ-Dataset-Doe)** (version 2: 265 cases over all three parameters, converged to 1e-6, with a mesh study and hold-out set). The cases here were stopped at 1e-3. Rerun to 1e-6, their velocity fields differ by typically 1 to 5 % (up to about 11 % for a few crucible cases), and the temperature cases at 1765 to 1785 K turn out not to be steady. This dataset stays online unchanged as version 1.

33 steady CFD cases of a Czochralski (CZ) silicon melt from three parameter
sweeps, collected in one place. I made them as reference data for surrogate
models and physics-informed neural networks (PINNs) in the Quasi project
(Aditya Seshaditya, Berlin).

The POD + RBF surrogate built on this data is in
[CZ_Surrogate_App](https://github.com/bertwinshine/CZ_Surrogate_App).

Case setup, CFD runs, mesh study and certification: Bertwin Kurisinkal Shine,
TU Chemnitz, 2025 to 2026.

## Setup

Geometry and boundary conditions follow:

> Huang, W., Cao, S., Liu, R., & Liu, D. (2025). Research on the thermal-fluid coupling in the growth process of Czochralski silicon single crystals based on an improved physics-informed neural network. *AIP Advances*, 15(10), 105202. https://doi.org/10.1063/5.0271778

All cases were run in ANSYS Fluent 2026 R1 as 2D axisymmetric steady cases of
the melt (radius 0.30 m, height 0.15 m). Every case is exported on the same
8181-node mesh, in the same node order, with the columns
`r, z, u_r, u_z, u_swirl, p, T`.

## Sweeps

One parameter changes at a time. The other two stay at the shared anchor:
crystal 8 rpm, crucible −3 rpm, hot wall 1745 K.

| Sweep | Folder | Range | Cases | Main result |
|---|---|---|---:|---|
| Crystal rotation | `CZ_Crystal_Sweep/` | 4 to 20 rpm | 11 | meridional flow reverses between 6 and 7 rpm |
| Hot wall temperature | `CZ_TempChange_Sweep/` | 1730 to 1785 K | 10 | flow stays steady up to 1785 K |
| Crucible rotation | `CZ_Crucible_Sweep/` | −1 to −10 rpm | 12 | velocity minimum near −4 rpm |

The anchor case appears in all three sweeps. The three copies are the same field
and differ only by export rounding (largest difference 5 × 10⁻⁶).

## Certification

Every case passed four checks before it was included:

1. Convergence: residuals below the set threshold with no drift.
2. Mesh independence: results stable against a refined mesh.
3. Conservation: mass and energy balance over the domain.
4. Boundary conditions: `certify_case.py` reads the swept value from the file
   name, works out every boundary condition from it (wall swirl, free-surface
   ramp, temperatures) and checks the exported field against them. It also
   checks that radial and swirl velocity are zero on the axis.

To rerun check 4 on a sweep:

```
cd CZ_Crystal_Sweep
python certify_case.py *.csv
```

## Structure

```
CZ_Crystal_Sweep/       11 cases + certify_case.py
CZ_Crucible_Sweep/      12 cases + certify_case.py
CZ_TempChange_Sweep/    10 cases + certify_case.py
```

Each sweep was first published as its own repository, with boundary-condition
tables, results, certification margins and an archived release with a DOI:

| Sweep | Repository | DOI |
|---|---|---|
| Crystal rotation | [CZ_Crystal_Sweep](https://github.com/bertwinshine/CZ_Crystal_Sweep) | [10.5281/zenodo.21955299](https://doi.org/10.5281/zenodo.21955299) |
| Hot wall temperature | [CZ_Study_TempChange](https://github.com/bertwinshine/CZ_Study_TempChange) | [10.5281/zenodo.21955315](https://doi.org/10.5281/zenodo.21955315) |
| Crucible rotation | [CZ_study_Crucible_Sweep](https://github.com/bertwinshine/CZ_study_Crucible_Sweep) | [10.5281/zenodo.21955323](https://doi.org/10.5281/zenodo.21955323) |

Not included here from the temperature repository:

- 1750 K, which was run on an earlier mesh and doesn't stack with the others.
- 1790, 1795 and 1800 K. The melt is weakly unsteady there, so those cases were
  run transient and exported as time averages. They are a different kind of data.

## Use

All 33 cases share the same nodes in the same order, so they can be stacked
directly as a POD snapshot matrix or used as reference fields for PINN
residuals.

If you use the data, please cite Huang et al. (2025), which the setup is based on.

## Limitations

- ANSYS Student licence (mesh capped at about 512k cells). The 8181-node 2D mesh
  was checked for mesh independence within that limit.
- Laminar throughout. The Reynolds number of the melt flow is well below
  transition for this geometry.
- The free surface is a prescribed wall, not solved. No turbulence model,
  radiation, Marangoni convection or magnetic field.
- One geometry only (the Huang et al. baseline). No changes in crucible shape,
  melt height or material properties.
- Only one parameter moves at a time, so the data says nothing about how the
  three interact.
- The fixed anchor rotations were entered as 0.838 rad/s for the crystal
  (8.0023 rpm) and −0.3140 rad/s for the crucible (−2.9985 rpm), 0.03 % and
  0.05 % off nominal. All swept values are exact.
- The certify scripts read the swept value from the file name. Keep the naming
  pattern if you rename files.

## Credit

The ANSYS setup, sweep design and validation are my work. I wrote the Python
with help from an AI assistant and have checked and understand all of it.
