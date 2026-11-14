#!/usr/bin/env python3
"""
Check CZ crystal-rotation cases against their boundary conditions.

The crystal rpm is read from the file name and the boundary conditions are
worked out from it. The free-surface swirl ramp is also fitted back to its
inner end, because a ramp left at a neighbouring case's rpm still looks
physical when compared point by point.

Usage:
    python certify_case.py steady/CZ_crystal_08rpm_steady.csv
    python certify_case.py steady/*.csv

Expected schema: r, z, u_r, u_z, u_swirl, p, T   (8181 rows)

Bertwin Kurisinkal Shine, TU Chemnitz. Written with help from an AI assistant.
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
OM_CRU  = -0.3140     # crucible rotation (rad/s), fixed in this sweep
T_CRY   = 1685.0      # crystal face temperature (K)
T_HOT   = 1745.0      # crucible wall temperature (K)
T_BOT_0 = 1700.0      # crucible bottom temperature at the axis (K)
RPM2RAD = 2 * np.pi / 60

TOL_REL   = 1e-3      # relative tolerance on prescribed wall swirl (0.1 %)
TOL_FIT   = 5e-3      # relative tolerance on the fitted free-surface inner end
TOL_T     = 0.05      # boundary temperature tolerance (K), away from corners
TOL_AXIS  = 1.0       # max |u_r| on axis, as % of field max |u_r|
N_NODES   = 8181


def rpm_from_name(path):
    """CZ_crystal_08rpm_steady.csv -> 8.0"""
    m = re.search(r"_(\d+(?:p\d+)?)rpm", os.path.basename(path))
    if not m:
        raise ValueError("cannot read rpm from filename: %s" % path)
    return float(m.group(1).replace("p", "."))


def certify(path):
    rpm = rpm_from_name(path)
    om = rpm * RPM2RAD
    v_inner = om * R_CRY                       # swirl at the crystal edge
    v_outer = OM_CRU * R_CRU                   # swirl at the crucible wall
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

    # --- swept boundary: the crystal face ----------------------------------
    # node export interpolates wall values, so use a relative tolerance
    ci = crys.r < R_CRY - 1e-6
    e = np.abs(crys.u_swirl - om * crys.r)[ci].max()
    add("crystal_face swirl", e / v_inner < 5e-3,
        "max err %.2e  (%.4f %% of edge speed %.6f)" % (e, 100 * e / v_inner, v_inner))

    # --- swept boundary: the free-surface ramp -----------------------------
    # point by point, then fitted back to the inner end (the fit catches a
    # ramp left at the wrong rpm)
    inner = free.r > R_CRY + 1e-6
    ramp = v_inner + ((free.r - R_CRY) / (R_CRU - R_CRY)) * (v_outer - v_inner)
    e = np.abs(free.u_swirl - ramp)[inner].max()
    add("free_surface swirl ramp", e < 1e-3, "max err %.2e" % e)

    mid = (free.r > R_CRY + 1e-6) & (free.r < R_CRU - 1e-6)
    slope, icept = np.polyfit(free.r[mid], free.u_swirl[mid], 1)
    fit_in = slope * R_CRY + icept
    rel = abs(fit_in - v_inner) / abs(v_inner)
    add("free_surface inner end", rel < TOL_FIT,
        "fitted %.6f  target %.6f  (%.3f %%)" % (fit_in, v_inner, 100 * rel))

    # --- fixed boundary: the crucible --------------------------------------
    e = np.abs(wall.u_swirl - v_outer).max()
    add("crucible_wall swirl", e / abs(v_outer) < TOL_REL,
        "max err %.2e = %.4f %%  (target %+.6f)" % (e, 100 * e / abs(v_outer), v_outer))

    e = np.abs(bot.u_swirl - OM_CRU * bot.r).max()
    add("crucible_bottom swirl", e / abs(v_outer) < TOL_REL,
        "max err %.2e = %.4f %%  (omega %+.6f)" % (e, 100 * e / abs(v_outer), OM_CRU))

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
    # u_r must be zero on the axis; also catches swapped u_r / u_z columns
    ratio = 100 * np.abs(axis.u_r).max() / np.abs(d.u_r).max()
    add("axis symmetry (u_r -> 0)", ratio < TOL_AXIS,
        "%.3f %% of field max" % ratio)

    ratio_s = 100 * np.abs(axis.u_swirl).max() / np.abs(d.u_swirl).max()
    add("axis swirl -> 0", ratio_s < TOL_AXIS, "%.3f %% of field max" % ratio_s)

    add("no NaN / inf", np.isfinite(d.values).all(), "")

    # --- report -------------------------------------------------------------
    passed = all(c[1] for c in checks)
    print("=" * 72)
    print("%s   (crystal %+g rpm)" % (os.path.basename(path), rpm))
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
