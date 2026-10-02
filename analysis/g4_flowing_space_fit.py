#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Gate G4: fit the flowing-space model (Formulation A) to SPARC.

Models (both 2 free parameters, omega and R_t, fit to V_obs):
  RT-additive : V = Vb + u
  Flow-A      : V = u - S R/2 + sqrt(S^2 R^2/4 + Vb^2)   (derived)
with u = omega R / (1 + R/R_t) and S = u/R - du/dR.

NFW and MOND-free BICs come from the 2026b tournament table. The RT-additive
refit is checked against that table as a pipeline sanity test.
Pre-registered in preregistration/PREREG.md.
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import models  # noqa: E402
from src.data import load_mass_models, load_tournament, UPSILON_DISK, UPSILON_BULGE  # noqa

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")


def vb2_signed(g):
    return (np.abs(g.Vgas) * g.Vgas + UPSILON_DISK * np.abs(g.Vdisk) * g.Vdisk
            + UPSILON_BULGE * np.abs(g.Vbul) * g.Vbul).values


def make_models(R, vb2):
    vb = np.sign(vb2) * np.sqrt(np.abs(vb2))

    def additive(_, om, Rt):
        return vb + models.rt_flow(R, om, Rt)

    def flow_a(_, om, Rt):
        u = models.rt_flow(R, om, Rt)
        S = models.rt_flow_shear(R, om, Rt)
        return u - S * R / 2 + np.sqrt(np.maximum((S * R) ** 2 / 4 + vb2, 0))

    return {"RT_additive": additive, "Flow_A": flow_a}


def fit_one(f, R, v, e):
    rmax = R.max()
    best = None
    for p0 in [(10, rmax / 2), (30, 2.0), (5, rmax * 2), (60, 0.5)]:
        try:
            p, cov = curve_fit(f, R, v, p0=p0, sigma=e, absolute_sigma=True,
                               bounds=([0, 0.1], [200, 5 * rmax]), maxfev=20000)
        except (RuntimeError, ValueError):
            continue
        chi2 = float(np.sum(((v - f(R, *p)) / e) ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, p, np.sqrt(np.diag(cov)))
    return best


def main():
    canon = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv",
                        index_col=0)
    tour = load_tournament()
    mm = load_mass_models()
    rows = []
    for gid, g in mm.groupby("ID"):
        if gid not in canon.index or canon.loc[gid, "Fit_Flag"] != "OK":
            continue
        R, v, e = g.R.values, g.Vobs.values, g.e_Vobs.values
        n = len(R)
        row = {"GalaxyID": gid, "N": n, "clean": bool(canon.loc[gid, "clean"])}
        for name, f in make_models(R, vb2_signed(g)).items():
            b = fit_one(f, R, v, e)
            if b is None:
                row[f"{name}_BIC"] = np.nan
                continue
            chi2, p, pe = b
            row.update({f"{name}_omega": p[0], f"{name}_Rt": p[1],
                        f"{name}_chi2": chi2,
                        f"{name}_BIC": chi2 + 2 * np.log(n)})
        if gid in tour.index:
            for m in ["RT", "NFW", "MondFree"]:
                row[f"tour_{m}_BIC"] = tour.loc[gid, f"{m}_BIC"]
        rows.append(row)
    df = pd.DataFrame(rows).set_index("GalaxyID")
    df.to_csv(ROOT / "results" / "g4_fits.csv")

    # Sanity: our RT-additive refit vs the tournament's RT fit.
    d = (df["RT_additive_BIC"] - df["tour_RT_BIC"]).abs()
    print("=" * 72)
    print("Gate G4: flowing-space (Formulation A) fits")
    print("=" * 72)
    print(f"\nPipeline check: |BIC_refit - BIC_tournament| for RT-additive:"
          f" median {d.median():.2f}, 90th pct {d.quantile(0.9):.2f}")

    res = {}
    for label, s in [("all_ok", df), ("clean", df[df["clean"]])]:
        out = {"N": int(len(s))}
        comps = {
            "Flow_A - RT_additive": s["Flow_A_BIC"] - s["RT_additive_BIC"],
            "Flow_A - NFW": s["Flow_A_BIC"] - s["tour_NFW_BIC"],
            "Flow_A - MOND_free": s["Flow_A_BIC"] - s["tour_MondFree_BIC"],
            "RT_additive - NFW": s["RT_additive_BIC"] - s["tour_NFW_BIC"],
        }
        print(f"\n[{label}] N = {len(s)}   (negative dBIC favours the first model)")
        for k, x in comps.items():
            x = x.dropna()
            out[k] = dict(median=float(x.median()),
                          frac_first_better=float((x < -2).mean()),
                          frac_second_better=float((x > 2).mean()))
            print(f"  {k:22s}: median dBIC {x.median():+7.2f};"
                  f" first better {100*(x < -2).mean():4.1f}%,"
                  f" second better {100*(x > 2).mean():4.1f}%")
        bics = s[["Flow_A_BIC", "RT_additive_BIC", "tour_NFW_BIC",
                  "tour_MondFree_BIC"]].dropna()
        wins = bics.idxmin(axis=1).value_counts(normalize=True)
        out["win_fractions"] = wins.to_dict()
        print("  outright BIC winners: " + ", ".join(
            f"{k.replace('_BIC', '').replace('tour_', '')} {100*v:.0f}%"
            for k, v in wins.items()))
        med_nfw = out["Flow_A - NFW"]["median"]
        out["competitive_with_NFW"] = bool(med_nfw <= 2)
        print(f"  pre-registered: Flow_A competitive with NFW"
              f" (median dBIC <= 2)? {'YES' if med_nfw <= 2 else 'NO'}")
        res[label] = out

    (ROOT / "results" / "g4_flowing_space_fit.json").write_text(
        json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
