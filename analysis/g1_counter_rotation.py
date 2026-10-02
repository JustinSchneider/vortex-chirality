#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Gate G1: counter-rotation chirality test.

For two counter-rotating disks at radius R, H_vort predicts drift-corrected
speeds differing by R*zeta (derivations/01_force_laws.py, section C):
    R*zeta_hat = V1 - V2 - k (s2^2 - s1^2) / (V1 + V2)
Expected |R*zeta| comes from clean SPARC galaxies of similar V_flat,
using the vorticity of their fitted RT flow at the same R.

Systems are listed in data/counter_rotators.csv with sources. Systems marked
blind = False were viewed before pre-registration and do not count toward
the overall kill criterion (preregistration/PREREG.md).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import models  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RNG = np.random.default_rng(20261002)
N_MC = 100_000
K_VALUES = (1.0, 2.0, 3.0)
DV_FLAT = 30.0   # km/s window for comparable SPARC galaxies


def expected_distribution(canon, V_sys, R):
    c = canon[canon["clean"] & (canon["Vflat"] > 0)
              & ((canon["Vflat"] - V_sys).abs() < DV_FLAT)]
    rz = R * models.rt_flow_vorticity(R, c["omega"].values, c["Rt"].values)
    return rz, len(c)


def main():
    canon = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv",
                        index_col=0)
    systems = pd.read_csv(ROOT / "data" / "counter_rotators.csv", comment="#")
    print("=" * 72)
    print("Gate G1: counter-rotation chirality")
    print("=" * 72)
    results = []
    for _, s in systems.iterrows():
        sini = np.sin(np.radians(s["inc_deg"]))
        V1 = RNG.normal(s["V1"], s["e_V1"], N_MC) / sini
        V2 = RNG.normal(s["V2"], s["e_V2"], N_MC) / sini
        s1 = np.abs(RNG.normal(s["sig1"], s["e_sig1"], N_MC))
        s2 = np.abs(RNG.normal(s["sig2"], s["e_sig2"], N_MC))
        V_sys = 0.5 * (s["V1"] + s["V2"]) / sini

        rz_exp, n_ref = expected_distribution(canon, V_sys, s["R_kpc"])
        p16 = float(np.percentile(rz_exp, 16))
        med = float(np.median(rz_exp))

        print(f"\n{s['system']} ({s['source']}), R = {s['R_kpc']} kpc,"
              f" i = {s['inc_deg']}°, blind = {s['blind']}")
        print(f"  expected R*zeta from {n_ref} SPARC galaxies with"
              f" V_flat = {V_sys:.0f} ± {DV_FLAT:.0f}:"
              f" median {med:.0f}, 16th pct {p16:.0f} km/s")
        upper_all = []
        per_k = {}
        kind = s.get("kind", "same_radius")
        for k in K_VALUES:
            if kind == "reversal":
                # Both sides are cold gas: R*zeta = |V_out| - |V_in| at R_rev.
                rz = V1 - V2
            else:
                k1 = 0.0 if s.get("cold1", False) else k
                k2 = 0.0 if s.get("cold2", False) else k
                rz = models.chirality_estimator(V1, V2, s1, s2, k1, k2)
            lo, mid, hi = np.percentile(rz, [2.3, 50, 97.7])
            upper = float(np.percentile(np.abs(rz), 97.7))
            upper_all.append(upper)
            per_k[k] = dict(median=float(mid), lo2s=float(lo), hi2s=float(hi),
                            abs_upper_2s=upper)
            print(f"  k = {k:.0f}: R*zeta_hat = {mid:+6.1f}"
                  f"  (2σ: {lo:+6.1f} to {hi:+6.1f});  |.| < {upper:5.1f}")
            if kind == "reversal":
                break
        bound = max(upper_all)
        disfavoured = bound < p16
        print(f"  max 2σ bound |R*zeta| < {bound:.0f} km/s vs 16th pct"
              f" {p16:.0f} -> {'DISFAVOURS H_vort' if disfavoured else 'not decisive'}")
        results.append(dict(system=s["system"], source=s["source"],
                            blind=bool(s["blind"]), R_kpc=float(s["R_kpc"]),
                            inc=float(s["inc_deg"]), per_k=per_k,
                            bound_2s=bound, expected_median=med,
                            expected_p16=p16, n_ref=n_ref,
                            disfavoured=bool(disfavoured)))

    blind = [r for r in results if r["blind"]]
    n_dis = sum(r["disfavoured"] for r in blind)
    print("\n" + "-" * 72)
    print(f"Blind systems: {len(blind)}; disfavouring H_vort: {n_dis}")
    if len(blind) == 0:
        print("Overall kill criterion not evaluable: no blind systems yet.")
    (ROOT / "results" / "g1_counter_rotation.json").write_text(
        json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
