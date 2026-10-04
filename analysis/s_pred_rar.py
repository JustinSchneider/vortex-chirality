#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""POST-HOC robustness check (not registered): the predicted chirality
amplitude normalised to the Radial Acceleration Relation instead of the
registered SPARC flow-profile fits.

If the velocity-dependent term carries the whole mass discrepancy, then in
the flowing-space picture with locally solid-body flow V_obs = V_bar + u, so
u(R) = V_RAR(R) - V_bar(R), where V_RAR follows from the RAR
(McGaugh, Lelli & Schombert 2016):
    g_obs = g_bar / (1 - exp(-sqrt(g_bar / g_dagger))),  g_dagger = 1.2e-10 m/s^2.
The asymmetry is R*zeta = d(R u)/dR. It is evaluated for every SPARC galaxy
(Q < 3, i >= 30) at the MaNGA counter-rotator's R_e, using SPARC galaxies with
|V_flat - V_c,*| < 30 km/s, exactly as in the registered S_pred.

SUPERSEDED (PREREG.md, 2026-10-04): the additive relation V_obs = V_bar + u
holds only for solid-body flow. For a general flow v_pro = V_bar + u - S R/2
+ O(S^2), so the first-order-consistent prediction is Delta V = 2 (V_RAR -
V_bar) (derivations/01_force_laws.py [A6]; analysis/agm_comparison.py "n2").
This script and pred_rz are kept unchanged so that the registered C1 numbers
(P = 44.9 km/s) remain reproducible; the paper no longer uses this formula.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.data import load_mass_models, load_sparc_table, v_bary  # noqa: E402

G_DAGGER = 3703.0       # km^2 s^-2 kpc^-1  (1.2e-10 m/s^2)


def rar_profiles():
    sp = load_sparc_table()
    sp = sp[(sp.Q < 3) & (sp.Inc >= 30) & (sp.Vflat > 0)]
    mm = load_mass_models()
    prof = {}
    for gid, g in mm.groupby("ID"):
        if gid not in sp.index:
            continue
        R = g.R.values
        vb = np.clip(v_bary(g.Vgas.values, g.Vdisk.values, g.Vbul.values), 0, None)
        gbar = vb ** 2 / R
        gobs = gbar / (1 - np.exp(-np.sqrt(np.clip(gbar, 1e-6, None) / G_DAGGER)))
        u = np.sqrt(gobs * R) - vb
        Ru = R * u
        rz = np.gradient(Ru, R)                      # d(R u)/dR = R * zeta
        prof[gid] = (R, rz, sp.loc[gid, "Vflat"])
    return prof


def main():
    prof = rar_profiles()
    t = pd.read_csv(ROOT / "results" / "g1_manga_sample.csv", index_col=0)
    cr = t[t.group == "CR"]
    vals, used = [], 0
    for _, g in cr.iterrows():
        rz = [np.interp(g.Re_kpc, R, z) for R, z, vf in prof.values()
              if abs(vf - g.Vc_star) < 30 and R.min() <= g.Re_kpc <= R.max()]
        if rz:
            vals.append(np.mean(rz))
            used += 1
    S_rar = float(np.mean(vals))
    reg = json.load(open(ROOT / "results" / "g1_manga.json"))
    print(f"SPARC galaxies with RAR profiles: {len(prof)}")
    print(f"S_pred (registered, flow-profile template) = {reg['S_pred']:.1f} km/s")
    print(f"S_pred (post-hoc, RAR normalisation)       = {S_rar:.1f} km/s  (from {used} CR galaxies)")
    p = reg["radii"]["1.0"]
    for k in ("D_lo", "D_mid", "D_hi"):
        S, s = p[k]["S"], p[k]["sS"]
        print(f"  {k}: S = {S:+.1f} ± {s:.1f}; S_RAR excluded at {(S_rar - S) / s:.1f} sigma;"
              f" 2-sigma upper bound = {100 * (S + 2 * s) / S_rar:.0f}% of S_RAR")
    (ROOT / "results" / "s_pred_rar.json").write_text(
        json.dumps({"S_pred_RAR": S_rar, "n_cr": used, "n_sparc": len(prof)}, indent=2))


if __name__ == "__main__":
    main()
