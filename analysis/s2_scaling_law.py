#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""S2a/S2b: which form can the scaling law R_t/R_d = F(Pi2, Pi3) take?

S2a: a theory with no disk length predicts R_t ∝ M_bar^0.5 and
     omega * V_sat / a0 = const. Rejected if the slope differs from 0.5 by >3 sigma.
S2b: F depends on Pi2 = G M_bar / (a0 R_d^2). Supported if permutation p < 0.01.

Pre-registered in preregistration/PREREG.md. Requires step0_foundation.py.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RNG = np.random.default_rng(20261002)
N_PERM = 10_000


def perm_spearman(x, y):
    rho = stats.spearmanr(x, y).statistic
    y = np.asarray(y)
    null = np.array([stats.spearmanr(x, RNG.permutation(y)).statistic
                     for _ in range(N_PERM)])
    return rho, (np.sum(np.abs(null) >= abs(rho)) + 1) / (N_PERM + 1)


def main():
    df = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv",
                     index_col=0)
    c = df[df["clean"]].copy()
    res = {"N": int(len(c))}
    print("=" * 72)
    print(f"S2: scaling-law form (clean sample, N = {len(c)})")
    print("=" * 72)

    # S2a
    r = stats.linregress(np.log10(c["M_bar"]), np.log10(c["Rt"]))
    z = (r.slope - 0.5) / r.stderr
    res["S2a"] = dict(slope=r.slope, se=r.stderr, z_vs_half=z,
                      rejected=bool(abs(z) > 3))
    q = np.percentile(c["omegaVsat_over_a0"], [16, 50, 84])
    res["S2a"]["omegaVsat_over_a0_16_50_84"] = q.tolist()
    rho_w, p_w = perm_spearman(np.log10(c["M_bar"]),
                               np.log10(c["omegaVsat_over_a0"]))
    res["S2a"]["omegaVsat_vs_Mbar"] = dict(rho=rho_w, p_perm=p_w)
    print(f"\n[S2a] log Rt = ({r.slope:.3f} ± {r.stderr:.3f}) log M_bar + c")
    print(f"      pure-a0 theory requires 0.5: z = {z:.1f}"
          f" -> {'REJECTED' if abs(z) > 3 else 'not rejected'}")
    print(f"      omega*V_sat/a0: 16/50/84% = {q[0]:.2f} / {q[1]:.2f} / {q[2]:.2f}"
          f"; trend with M_bar rho = {rho_w:.2f} (p = {p_w:.1e})")

    # S2b
    rho, p = perm_spearman(np.log10(c["Pi2"]), np.log10(c["Rt_over_Rd"]))
    fit = stats.linregress(np.log10(c["Pi2"]), np.log10(c["Rt_over_Rd"]))
    res["S2b"] = dict(rho=rho, p_perm=p, slope=fit.slope, se=fit.stderr,
                      r2=fit.rvalue ** 2, supported=bool(p < 0.01))
    print(f"\n[S2b] Spearman rho(Rt/Rd, Pi2) = {rho:.3f}, permutation p = {p:.1e}"
          f" -> {'SUPPORTED' if p < 0.01 else 'not supported'}")
    print(f"      power law: Rt/Rd ∝ Pi2^({fit.slope:.3f} ± {fit.stderr:.3f}),"
          f" R² = {fit.rvalue**2:.3f}")
    print(f"      Pi2 range: {c['Pi2'].min():.2f} - {c['Pi2'].max():.1f}"
          f" (median {c['Pi2'].median():.2f})")

    # Exploratory (not pre-registered): B/T as a stand-in for Pi3, partial
    # on Pi2, among galaxies with bulges. Labelled exploratory in output.
    from_rank = lambda s: stats.rankdata(s)
    x, y, zc = (from_rank(c["BT"]), from_rank(c["Rt_over_Rd"]),
                from_rank(c["Pi2"]))
    def resid(a, b):
        return a - np.polyval(np.polyfit(b, a, 1), b)
    pr = stats.pearsonr(resid(x, zc), resid(y, zc))
    res["exploratory_BT_partial_Pi2"] = dict(rho=pr.statistic, p=pr.pvalue,
                                             n_bulged=int((c["BT"] > 0).sum()))
    print(f"\n[exploratory] partial rho(Rt/Rd, B/T | Pi2) = {pr.statistic:.3f}"
          f" (p = {pr.pvalue:.3f}); bulged galaxies in clean sample:"
          f" {(c['BT'] > 0).sum()}")

    (ROOT / "results" / "s2_scaling_law.json").write_text(
        json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
