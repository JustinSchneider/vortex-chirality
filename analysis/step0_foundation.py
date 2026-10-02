#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 0: rebuild the empirical foundation from the existing RT fits.

Builds data/processed/sparc_canonical.csv (one row per galaxy, with
cleaning flags and the Pi groups) and recomputes the headline RT
statistics with correct uncertainties. See preregistration/PREREG.md.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.data import (load_bulges, load_mass_models, load_rt_fits,  # noqa
                      load_sparc_table, UPSILON_DISK, UPSILON_BULGE, HELIUM)

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

G_KPC = 4.30091e-6          # kpc (km/s)^2 / Msun
A0 = 3703.0                 # km^2 s^-2 kpc^-1 (1.2e-10 m/s^2)
RNG = np.random.default_rng(20261002)
N_PERM = 10_000


def build_table():
    sp = load_sparc_table()
    fits = load_rt_fits()
    lbul = load_bulges().reindex(sp.index).fillna(0.0)
    mm = load_mass_models()
    rmax = mm.groupby("ID")["R"].max().rename("R_max")

    df = sp.join(fits, how="inner").join(rmax)
    df["Lbul"] = lbul
    df["Ldisk"] = df["L36"] - df["Lbul"]
    df["M_bar"] = 1e9 * (HELIUM * df["MHI"] + UPSILON_DISK * df["Ldisk"]
                         + UPSILON_BULGE * df["Lbul"])
    df["BT"] = df["Lbul"] / df["L36"]
    df["omega"] = df["Tapered_omega"]
    df["Rt"] = df["Tapered_Rt"]
    df["Rt_err"] = df["Tapered_Rt_err"]
    df["V_sat"] = df["omega"] * df["Rt"]
    df["Rt_over_Rd"] = df["Rt"] / df["Rdisk"]
    df["Pi2"] = G_KPC * df["M_bar"] / (A0 * df["Rdisk"] ** 2)
    df["omegaVsat_over_a0"] = df["omega"] * df["V_sat"] / A0

    df["ok_quality"] = df["Q"] < 3
    df["ok_inc"] = df["Inc"] >= 30
    df["ok_fit"] = df["Fit_Flag"] == "OK"
    df["ok_ident"] = (df["Rt_err"] < df["Rt"]) & (df["Rt"] < df["R_max"])
    df["clean"] = df[["ok_quality", "ok_inc", "ok_fit", "ok_ident"]].all(axis=1)
    df["taper_pref"] = df["ok_fit"] & (df["Delta_BIC"] > 2)
    return df


def ols_loglog(x, y):
    r = stats.linregress(np.log10(x), np.log10(y))
    return dict(slope=r.slope, slope_se=r.stderr, intercept=r.intercept,
                r2=r.rvalue ** 2, p=r.pvalue, n=len(x))


def perm_p_spearman(x, y, n=N_PERM):
    rho = stats.spearmanr(x, y).statistic
    y = np.asarray(y)
    null = np.array([stats.spearmanr(x, RNG.permutation(y)).statistic
                     for _ in range(n)])
    return rho, (np.sum(np.abs(null) >= abs(rho)) + 1) / (n + 1)


def main():
    df = build_table()
    out = ROOT / "data" / "processed" / "sparc_canonical.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out)

    res = {}
    print("=" * 72)
    print("Step 0: foundation")
    print("=" * 72)
    print(f"\nGalaxies: {len(df)}")
    for flag in ["ok_quality", "ok_inc", "ok_fit", "ok_ident"]:
        print(f"  fail {flag:11s}: {(~df[flag]).sum():3d}")
    print(f"  clean sample : {df['clean'].sum()}")
    print(f"  taper-preferred (dBIC>2, fit OK): {df['taper_pref'].sum()}")
    print(f"  low-inclination galaxies inside the paper's 171:"
          f" {(df['ok_fit'] & ~df['ok_inc']).sum()}")

    samples = {
        "paper_171": df[df["ok_fit"]],
        "taper_pref": df[df["taper_pref"]],
        "clean": df[df["clean"]],
    }
    for name, s in samples.items():
        reg = ols_loglog(s["Rdisk"], s["Rt"])
        rho, p_perm = perm_p_spearman(np.log10(s["Rdisk"]), np.log10(s["Rt"]))
        med = float(np.median(s["Rt_over_Rd"]))
        q16, q84 = np.percentile(s["Rt_over_Rd"], [16, 84])
        res[name] = dict(**reg, median_Rt_Rd=med, Rt_Rd_16=q16, Rt_Rd_84=q84,
                         spearman=rho, p_perm=p_perm)
        print(f"\n[{name}] N = {reg['n']}")
        print(f"  log Rt = ({reg['slope']:.3f} ± {reg['slope_se']:.3f}) log Rd"
              f" + {reg['intercept']:.3f};  R² = {reg['r2']:.3f}")
        print(f"  median Rt/Rd = {med:.2f}  (16-84%: {q16:.2f}-{q84:.2f})")
        print(f"  Spearman rho = {rho:.3f}, permutation p = {p_perm:.1e}")

    print("\nPaper 2026a quoted: slope 0.794 ± 0.063, R² = 0.135,"
          " median Rt/Rd = 2.42 (N = 171)")
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "step0_foundation.json").write_text(
        json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
