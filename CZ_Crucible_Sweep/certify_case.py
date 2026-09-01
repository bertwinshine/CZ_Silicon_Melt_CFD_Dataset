#!/usr/bin/env python3
"""
certify_case.py
===============

Checks one exported CZ crucible-sweep case against the boundary conditions it
should have been run with, plus the physical sanity tests.

The crucible rotation is read from the filename and every boundary condition is
re-derived from it, so a case whose setup does not match its label fails here
rather than entering the dataset unnoticed.

Usage:
    python certify_case.py steady/CZ_crucible_m5rpm_steady.csv
    python certify_case.py steady/*.csv

Expected schema: r, z, u_r, u_z, u_swirl, p, T   (8181 rows)

Author: Bertwin Kurisinkal Shine
"""

import sys
import os
import re
import numpy as np
import pandas as pd

# ---- geometry and fixed conditions ----------------------------------------
R_CRU   = 0.30        # crucible radius (m)
R_CRY   = 0.18        # crystal radius (m)
Z_TOP   = 0.15        # melt height (m)
OM_CRY  = 0.837758    # crystal rotation, +8 rpm (rad/s) -- fixed in this sweep
T_CRY   = 1685.0      # crystal face temperature (K)
T_HOT   = 1745.0      # crucible wall temperature (K)
T_BOT_0 = 1700.0      # crucible bottom temperature at the axis (K)
RPM2RAD = 2 * np.pi / 60

TOL_BC    = 1e-5      # boundary velocity tolerance (m/s)
TOL_REL   = 1e-3      # relative tolerance on prescribed wall swirl (0.1 %)
TOL_T     = 0.05      # boundary temperature tolerance (K), away from corners
TOL_AXIS  = 1.0       # max |u_r| on axis, as % of field max |u_r|
N_NODES   = 8181


def rpm_from_name(path):
    """CZ_crucible_m5p5rpm_steady.csv -> -5.5"""
    m = re.search(r"m(\d+(?:p\d+)?)rpm", os.path.basename(path))
    if not m:
        raise ValueError("cannot read rpm from filename: %s" % path)
    return -float(m.group(1).replace("p", "."))


def certify(path):
    rpm = rpm_from_name(path)
    om = rpm * RPM2RAD
    v_outer = om * R_CRU                       # swirl at the crucible wall
    v_inner = OM_CRY * R_CRY                   # swirl at the crystal edge
    d = pd.read_csv(path)

    checks = []

    def add(name, ok, detail):
        checks.append((name, ok, detail))

    # --- schema -------------------------------------------------------------
    cols = ["r", "z", "u_r", "u_z", "u_swirl", "p", "T"]
    add("schema", list(d.columns) == cols, ", ".join(d.columns))
    add("node count", len(d) == N_NODES, "%d rows" % len(d))

    # --- zone masks ---------------------------------------------------------
    top = d[np.abs(d.z - Z_TOP) < 1e-9]
    free = top[top.r >= R_CRY - 1e-9]
    crys = top[top.r <= R_CRY + 1e-9]
    wall = d[(np.abs(d.r - R_CRU) < 1e-9) & (d.z > 1e-6) & (d.z < Z_TOP - 1e-6)]
    bot = d[(d.z < 1e-12) & (d.r > 1e-6) & (d.r < R_CRU - 1e-6)]
    axis = d[d.r < 1e-9]

    # --- swirl boundary conditions -----------------------------------------
    # Every prescribed swirl condition is re-derived from the case rpm, so a
    # boundary that does not follow the swept parameter is caught here. The
    # corner node at r = R_CRY is shared with the crystal face and is blended
    # by the solver, so it is excluded from the ramp checks.
    e = np.abs(wall.u_swirl - v_outer).max()
    add("crucible_wall swirl", e / abs(v_outer) < TOL_REL,
        "max err %.2e = %.4f %%  (target %+.6f)" % (e, 100 * e / abs(v_outer), v_outer))

    e = np.abs(bot.u_swirl - om * bot.r).max()
    add("crucible_bottom swirl", e / abs(v_outer) < TOL_REL,
        "max err %.2e = %.4f %%  (omega %+.6f)" % (e, 100 * e / abs(v_outer), om))

    inner = free.r > R_CRY + 1e-6
    ramp = v_inner + ((free.r - R_CRY) / (R_CRU - R_CRY)) * (v_outer - v_inner)
    e = np.abs(free.u_swirl - ramp)[inner].max()
    add("free_surface swirl ramp", e < 1e-3,
        "max err %.2e  (%.6f -> %.6f)" % (e, v_inner, v_outer))

    # Node-based export interpolates the wall value, so this one is judged
    # relative to the crystal edge speed rather than on an absolute tolerance.
    ci = crys.r < R_CRY - 1e-6
    e = np.abs(crys.u_swirl - OM_CRY * crys.r)[ci].max()
    add("crystal_face swirl", e / v_inner < 5e-3,
        "max err %.2e  (%.3f %% of edge speed)" % (e, 100 * e / v_inner))

    # --- thermal boundary conditions ---------------------------------------
    ramp = T_CRY + ((free.r - R_CRY) / (R_CRU - R_CRY)) * (T_HOT - T_CRY)
    e = np.abs(free["T"] - ramp)[inner].max()
    add("free_surface T ramp", e < TOL_T, "max err %.3f K" % e)

    ramp = T_BOT_0 + (bot.r / R_CRU) * (T_HOT - T_BOT_0)
    e = np.abs(bot["T"] - ramp).max()
    add("crucible_bottom T ramp", e < 0.2, "max err %.3f K" % e)

    ok = (np.abs(d["T"].min() - T_CRY) < TOL_T) and (np.abs(d["T"].max() - T_HOT) < TOL_T)
    add("T within bounds", ok, "%.2f - %.2f K" % (d["T"].min(), d["T"].max()))

    # --- physical sanity ----------------------------------------------------
    # On the symmetry axis the radial velocity must vanish. This is a
    # consequence of axisymmetry, not a numerical tolerance, and it also
    # catches a u_r / u_z column transposition at export time.
    ratio = 100 * np.abs(axis.u_r).max() / np.abs(d.u_r).max()
    add("axis symmetry (u_r -> 0)", ratio < TOL_AXIS,
        "%.3f %% of field max" % ratio)

    ratio_s = 100 * np.abs(axis.u_swirl).max() / np.abs(d.u_swirl).max()
    add("axis swirl -> 0", ratio_s < TOL_AXIS, "%.3f %% of field max" % ratio_s)

    add("no NaN / inf", np.isfinite(d.values).all(), "")

    # --- report -------------------------------------------------------------
    passed = all(c[1] for c in checks)
    print("=" * 72)
    print("%s   (crucible %+g rpm)" % (os.path.basename(path), rpm))
    print("=" * 72)
    for name, ok, detail in checks:
        print("  [%s] %-26s %s" % ("PASS" if ok else "FAIL", name, detail))
    print("  mean |u| = %.6f   p range = %.4f   T mean = %.4f"
          % (np.sqrt(d.u_r**2 + d.u_z**2).mean(),
             d.p.max() - d.p.min(), d["T"].mean()))
    print("  RESULT: %s\n" % ("CERTIFIED" if passed else "REJECTED"))
    return passed


if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        print(__doc__)
        sys.exit(1)
    results = [certify(f) for f in files]
    n = len(results)
    print("%d of %d certified." % (sum(results), n))
    sys.exit(0 if all(results) else 1)
