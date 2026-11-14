#!/usr/bin/env python3
"""
Check CZ temperature-sweep cases against their boundary conditions.

The hot wall temperature is read from the file name and the thermal boundary
conditions are worked out from it.

Usage:
    python certify_case.py "steady/cz_(temp 1760).csv"
    python certify_case.py steady/*.csv

Filename must contain the hot-wall temperature as a 4-digit number, e.g.
    cz_(temp 1760).csv  ->  T_hot = 1760 K
    cz_(temp 1750).csv  ->  T_hot = 1750 K

Expected schema: r, z, u_r, u_z, u_swirl, p, T   (8181 rows)

Rotation is fixed in this sweep (crystal +8 rpm, crucible -3 rpm) and the
swirl is checked against it.

Bertwin Kurisinkal Shine, TU Chemnitz. Written with help from an AI assistant.
"""

import sys
import os
import re
import numpy as np
import pandas as pd

# ---- geometry -------------------------------------------------------------
R_CRU   = 0.30        # crucible radius (m)
R_CRY   = 0.18        # crystal radius (m)
Z_TOP   = 0.15        # melt height (m)

# ---- fixed thermal references (Huang et al. 2025, Table I) ----------------
T_CRY   = 1685.0      # crystal face temperature (K), fixed
T_BOT_0 = 1700.0      # crucible bottom temperature at the axis (K)

# ---- fixed rotation for the temperature sweep -----------------------------
OM_CRY  =  0.837758   # crystal rotation, +8 rpm (rad/s)
OM_CRU  = -0.314000   # crucible rotation, -3 rpm (rad/s)

# ---- tolerances -----------------------------------------------------------
TOL_T      = 0.05     # boundary temperature tolerance (K), away from corners
TOL_T_RAMP = 0.2      # thermal ramp tolerance (K)
TOL_T_FACE = 0.15     # crystal face tolerance (K); nodes near the corner are
                      # blended by the solver and read up to ~0.1 K high
TOL_REL    = 1e-3     # relative tolerance on fixed wall swirl (0.1 %)
TOL_AXIS   = 1.0      # max |u_r| on axis, as % of field max |u_r|
N_NODES    = 8181


def temp_from_name(path):
    '''cz_(temp 1760).csv -> 1760.0 ; also matches cz_temp_1760_... '''
    m = re.search(r"(\d{4})", os.path.basename(path))
    if not m:
        raise ValueError("cannot read temperature from filename: %s" % path)
    return float(m.group(1))


def certify(path):
    t_hot = temp_from_name(path)
    v_outer = OM_CRU * R_CRU                 # swirl at the crucible wall (fixed)
    v_inner = OM_CRY * R_CRY                 # swirl at the crystal edge (fixed)
    d = pd.read_csv(path)

    checks = []
    def add(name, ok, detail):
        checks.append((name, ok, detail))

    # --- schema -------------------------------------------------------------
    cols = ["r", "z", "u_r", "u_z", "u_swirl", "p", "T"]
    add("schema", list(d.columns) == cols, ", ".join(d.columns))
    add("node count", len(d) == N_NODES, "%d rows" % len(d))

    # --- zone masks ---------------------------------------------------------
    top  = d[np.abs(d.z - Z_TOP) < 1e-9]
    free = top[top.r >= R_CRY - 1e-9]
    crys = top[top.r <= R_CRY + 1e-9]
    wall = d[(np.abs(d.r - R_CRU) < 1e-9) & (d.z > 1e-6) & (d.z < Z_TOP - 1e-6)]
    bot  = d[(d.z < 1e-12) & (d.r > 1e-6) & (d.r < R_CRU - 1e-6)]
    axis = d[d.r < 1e-9]
    inner = free.r > R_CRY + 1e-6

    # --- thermal boundary conditions (re-derived from the filename T_hot) ---
    # crucible side wall: fixed hot temperature
    e = np.abs(wall["T"] - t_hot).max()
    add("crucible_wall T (fixed)", e < TOL_T, "max err %.3f K  (target %.1f)" % (e, t_hot))

    # free surface: linear ramp T_cry -> T_hot over r = R_CRY .. R_CRU
    ramp = T_CRY + ((free.r - R_CRY) / (R_CRU - R_CRY)) * (t_hot - T_CRY)
    e = np.abs(free["T"] - ramp)[inner].max()
    add("free_surface T ramp", e < TOL_T_RAMP, "max err %.3f K  (%.1f -> %.1f)" % (e, T_CRY, t_hot))

    # crucible bottom: linear ramp T_bot_0 (axis) -> T_hot (rim) over r = 0 .. R_CRU
    ramp = T_BOT_0 + (bot.r / R_CRU) * (t_hot - T_BOT_0)
    e = np.abs(bot["T"] - ramp).max()
    add("crucible_bottom T ramp", e < TOL_T_RAMP, "max err %.3f K  (%.1f -> %.1f)" % (e, T_BOT_0, t_hot))

    # crystal face; corner node at r = R_CRY is blended by the solver, excluded
    ci = crys.r < R_CRY - 1e-6
    e = np.abs(crys["T"] - T_CRY)[ci].max()
    add("crystal_face T (fixed)", e < TOL_T_FACE,
        "max err %.3f K  (target %.1f, corner node excluded)" % (e, T_CRY))

    # global bounds: coldest = crystal, hottest = hot wall
    ok = (np.abs(d["T"].min() - T_CRY) < TOL_T) and (np.abs(d["T"].max() - t_hot) < TOL_T)
    add("T within bounds", ok, "%.2f - %.2f K  (expect %.0f - %.0f)"
        % (d["T"].min(), d["T"].max(), T_CRY, t_hot))

    # --- rotation guard (fixed across this sweep) ---------------------------
    e = np.abs(wall.u_swirl - v_outer).max()
    add("crucible_wall swirl (fixed)", e / abs(v_outer) < TOL_REL,
        "max err %.2e = %.4f %%  (target %+.6f)" % (e, 100 * e / abs(v_outer), v_outer))

    e = np.abs(bot.u_swirl - OM_CRU * bot.r).max()
    add("crucible_bottom swirl (fixed)", e / abs(v_outer) < TOL_REL,
        "max err %.2e = %.4f %%  (omega %+.6f)" % (e, 100 * e / abs(v_outer), OM_CRU))

    # --- physical sanity ----------------------------------------------------
    # u_r must be zero on the axis; also catches swapped u_r / u_z columns
    ratio = 100 * np.abs(axis.u_r).max() / np.abs(d.u_r).max()
    add("axis symmetry (u_r -> 0)", ratio < TOL_AXIS, "%.3f %% of field max" % ratio)

    ratio_s = 100 * np.abs(axis.u_swirl).max() / np.abs(d.u_swirl).max()
    add("axis swirl -> 0", ratio_s < TOL_AXIS, "%.3f %% of field max" % ratio_s)

    # geometry not swapped
    add("domain size (r,z)",
        abs(d.r.max() - R_CRU) < 0.02 and abs(d.z.max() - Z_TOP) < 0.02,
        "r max %.4f, z max %.4f" % (d.r.max(), d.z.max()))

    add("no NaN / inf", np.isfinite(d.values).all(), "")
    add("no duplicate nodes", not d.duplicated(subset=["r", "z"]).any(),
        "%d duplicate (r,z)" % int(d.duplicated(subset=["r", "z"]).sum()))

    # --- report -------------------------------------------------------------
    passed = all(c[1] for c in checks)
    print("=" * 72)
    print("%s   (hot wall %g K)" % (os.path.basename(path), t_hot))
    print("=" * 72)
    for name, ok, detail in checks:
        print("  [%s] %-28s %s" % ("PASS" if ok else "FAIL", name, detail))
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
    print("%d of %d certified." % (sum(results), len(results)))
    sys.exit(0 if all(results) else 1)
