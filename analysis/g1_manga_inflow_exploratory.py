#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""EXPLORATORY (data already unblinded): does the counter-rotator gas
slowdown S track radial gas flow?

Conventional explanation: accreted counter-rotating gas loses angular
momentum and flows inward, so it rotates below circular speed, and the
slowdown should scale with |V_R|. H_vort predicts a slowdown that persists
at |V_R| = 0.

For the 81 CR + 243 matched CO galaxies of the registered test:
  1. Spearman rho(D, |V_R|/V_c) within each group;
  2. Theil-Sen fit of D on |V_R|/V_c in CR -> slowdown extrapolated to no
     radial flow, S0 = median(D_CO) - intercept;
  3. S restricted to galaxies with |V_R| below the combined median.
All at 1 R_e and the three k values. Not a registered test.
"""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import manga_measure as mm  # noqa: E402
from analysis.g1_manga_test import galaxy_arrays  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
RNG = np.random.default_rng(20261004)


def vr_one(args):
    pifu, pa_s, pa_g, inc, re = args
    maps, x, y, _, _ = galaxy_arrays(pifu)
    r = {row["r_re"]: row for row in mm.measure(maps, x, y, pa_s, pa_g, inc, re)}
    return pifu, r[1.0]["VRg"], r[1.5]["VRg"]


def boot_S(a, b, n=10_000):
    S = np.median(b) - np.median(a)
    bs = [np.median(RNG.choice(b, len(b))) - np.median(RNG.choice(a, len(a)))
          for _ in range(n)]
    return S, np.std(bs, ddof=1)


def main():
    t = pd.read_csv(ROOT / "results" / "g1_manga_sample.csv", index_col=0)
    jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in t.iterrows()]
    with ProcessPoolExecutor(12) as ex:
        vr = {p: (a, b) for p, a, b in ex.map(vr_one, jobs, chunksize=4)}
    t["VR1"] = [vr[p][0] for p in t.index]
    t["fR"] = t["VR1"] / t["Vc_star"]
    cr, co = t[t.group == "CR"], t[t.group == "CO"]
    out = {}
    print("EXPLORATORY: gas radial flow vs counter-rotator slowdown (1 R_e)")
    print(f"|V_R| median: CR {cr.VR1.median():.1f}, CO {co.VR1.median():.1f} km/s;"
          f"  |V_R|/V_c median: CR {cr.fR.median():.2f}, CO {co.fR.median():.2f}")
    vr_cut = t["VR1"].median()
    for k in ("D_lo", "D_mid", "D_hi"):
        col = f"{k}_1.0"
        a = cr[[col, "fR"]].dropna()
        b = co[[col, "fR"]].dropna()
        rho_cr = stats.spearmanr(a.fR, a[col])
        rho_co = stats.spearmanr(b.fR, b[col])
        ts = stats.theilslopes(a[col], a.fR)
        # bootstrap the extrapolated no-flow slowdown
        s0b = []
        for _ in range(2000):
            ia = RNG.integers(0, len(a), len(a))
            ib = RNG.integers(0, len(b), len(b))
            tsb = stats.theilslopes(a[col].values[ia], a.fR.values[ia])
            s0b.append(np.median(b[col].values[ib]) - tsb.intercept)
        S0 = np.median(b[col]) - ts.intercept
        lo = t[t.VR1 <= vr_cut]
        Sl, sSl = boot_S(lo[lo.group == "CR"][col].dropna().values,
                         lo[lo.group == "CO"][col].dropna().values)
        n_lo = int((lo.group == "CR").sum())
        out[k] = dict(rho_cr=rho_cr.statistic, p_cr=rho_cr.pvalue,
                      rho_co=rho_co.statistic, p_co=rho_co.pvalue,
                      slope=ts.slope, S0=S0, sS0=float(np.std(s0b, ddof=1)),
                      S_lowVR=Sl, sS_lowVR=sSl, n_cr_lowVR=n_lo)
        print(f"  {k:6s}: rho(D,|VR|/Vc) CR {rho_cr.statistic:+.2f} (p={rho_cr.pvalue:.3f}),"
              f" CO {rho_co.statistic:+.2f} (p={rho_co.pvalue:.3f});"
              f"  no-flow S0 = {S0:+5.1f} ± {np.std(s0b, ddof=1):4.1f};"
              f"  S(|VR| < {vr_cut:.0f}) = {Sl:+5.1f} ± {sSl:4.1f} (N_CR={n_lo})")
    (ROOT / "results" / "g1_manga_inflow_exploratory.json").write_text(
        json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
