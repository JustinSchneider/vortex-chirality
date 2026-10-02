#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Realistic power of the SIGNED test (flow follows the stellar disk), using
the same real nuisance properties as g1_manga_power_real.py. No gas speeds
or D values are used.

S = median(D_CO) - median(D_CR); H_vort with flow along the stars predicts
S ~ mean R*zeta > 0. Verdict rules mirror the registered ones:
  falsified: S + 2 sigma < 0.25 S_pred;  supported: S > 3 sigma and
  |S - S_pred| < 2 sigma;  both required at k-1, k, k+1.
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import analysis.g1_manga_power_real as P  # noqa: E402
from analysis.g1_manga_test import load_catalogues, match  # noqa: E402

RNG = np.random.default_rng(20261003)
N_SIM, N_BOOT = 300, 300


def boot_S(a, b):
    S = np.median(b) - np.median(a)
    bs = [np.median(RNG.choice(b, len(b))) - np.median(RNG.choice(a, len(a)))
          for _ in range(N_BOOT)]
    return S, np.std(bs, ddof=1)


def sim(cr, co, f, sp, signal):
    def nuis(g, ff):
        Vc = g["Vc_star"].values
        gain = (g["sg"].fillna(0).values ** 2 - g["ss"].values ** 2) / (2 * Vc)
        n = (ff * RNG.normal(0, P.GAS_FRAC * Vc) + RNG.normal(0, P.SIG_K, len(g)) * gain
             + RNG.normal(0, np.sqrt(g["eVs"].fillna(5) ** 2 + g["eVg"].fillna(5) ** 2)))
        return n, gain
    ncr, gcr = nuis(cr, f)
    nco, gco = nuis(co, 1.0)
    sig = np.zeros(len(cr))
    if signal:   # flow along stars: counter-rotating gas runs slow by R*zeta
        sig = -np.array([RNG.choice(P.rzeta_draw(sp, g["Vc_star"], g["Re_kpc"]))
                         for _, g in cr.iterrows()])
    delta = RNG.uniform(-1, 1)
    return [boot_S(sig + ncr + (j - delta) * gcr, nco + (j - delta) * gco) for j in (-1, 0, 1)]


def main():
    d = load_catalogues()
    d = d[(d.e_pa_star <= 20) & (d.e_pa_gas <= 20) & d.dpa.notna() & ~d.merger]
    d["group"] = np.where(d.dpa > 150, "CR", np.where(d.dpa < 30, "CO", ""))
    d = d[d.group != ""].copy()
    d["Re_kpc"] = d.Re_arcsec * d.kpc_per_arcsec
    d["log_sig1re"] = np.log10(d.sig1re.clip(lower=1))
    t = d.join(pd.read_csv(P.NUIS, index_col=0), how="inner")
    t = t[np.isfinite(t.Vs) & np.isfinite(t.ss) & (t.Vc_star > 0)].dropna(
        subset=["logMstar", "ttype", "log_sig1re", "inc"])
    t["resid_ratio"] = t.gas_resid / t.med_eg
    sp = P.sparc_pool()
    print("Signed test, realistic power (no gas speeds or D used)")
    for vcut in (100, 80, 60, 40):
        s = t[t.Vs > vcut]
        cr, pool = s[s.group == "CR"], s[s.group == "CO"]
        pairs = match(cr, pool, ["logMstar", "ttype", "log_sig1re", "inc"])
        co = pool.loc[[p for v in pairs.values() for p in v]]
        f = float(np.nanmedian(cr.resid_ratio) / np.nanmedian(co.resid_ratio))
        S_pred = float(np.mean([np.mean(P.rzeta_draw(sp, g.Vc_star, g.Re_kpc)) for _, g in cr.iterrows()]))
        for flab, ff in (("measured f", f), ("f x 1.3", 1.3 * f)):
            c = np.zeros((2, 2))
            for _ in range(N_SIM):
                for si, signal in enumerate((True, False)):
                    r = sim(cr, co, ff, sp, signal)
                    c[si, 0] += all(S + 2 * e < 0.25 * S_pred for S, e in r)
                    c[si, 1] += all(S > 3 * e and abs(S - S_pred) < 2 * e for S, e in r)
            p = c / N_SIM
            print(f"V_*>{vcut:3d} N_CR={len(cr):3d} S_pred={S_pred:4.0f} [{flab:10s}]"
                  f" signal: SUPPORTED {p[0,1]:.2f} falsified {p[0,0]:.2f} |"
                  f" null: FALSIFIED {p[1,0]:.2f} false-support {p[1,1]:.2f}", flush=True)


if __name__ == "__main__":
    main()
