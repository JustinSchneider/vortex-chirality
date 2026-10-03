#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""POST HOC (not pre-registered): counter-rotating gas against the strong-
gravitomagnetic GR model of Astesiano, Ruggiero & Re (2026, arXiv:2606.18868).

The model's prediction comes from derivations/04_strong_gravitomagnetic_branches.py:
  psi = C0 r  (their solution):   Delta V = (V_obs^2 - V_bar^2) / V_obs   ("n = 1")
  psi = K r^2 (rigid dragging):   Delta V = 2 (V_obs - V_bar)             ("n = 2")
evaluated with the RAR (g_dagger = 1.2e-10 m/s^2) on SPARC mass models, matched
to each counter-rotator by |V_flat - V_c,*| < 30 km/s and taken at the
counter-rotator's own radius, exactly as for the registered prediction P
(analysis/s_pred_rar.py). The model fits the observed (stellar) rotation curve
with its co-rotating branch, so the counter-rotating gas is on the
counter-rotating branch: the signed slowdown s is the relevant observable.

Statistic: the C1 primary statistic (zero-twist Theil-Sen intercept of
s_i = median(D of its own 3 controls) - D_i), with P replaced by the model
prediction. Also: SAMI descriptive values, an outer-disc forecast, and the
number of counter-rotators an outer-disc test would need.
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.c1_power import intercept_and_sd, log10_bf  # noqa: E402
from analysis.c1_test import primary  # noqa: E402
from analysis.g1_manga_robustness import COLS, build  # noqa: E402
from analysis.g1_manga_test import gas_twist, match  # noqa: E402
from src.data import load_mass_models, load_sparc_table, v_bary  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
G_DAGGER = 3703.0
RNG = np.random.default_rng(20261003)


def model_profiles():
    sp = load_sparc_table()
    sp = sp[(sp.Q < 3) & (sp.Inc >= 30) & (sp.Vflat > 0)]
    prof = {}
    for gid, g in load_mass_models().groupby("ID"):
        if gid not in sp.index:
            continue
        R = g.R.values
        vb = np.clip(v_bary(g.Vgas.values, g.Vdisk.values, g.Vbul.values), 0, None)
        gb = vb ** 2 / R
        vo = np.sqrt(gb / (1 - np.exp(-np.sqrt(np.clip(gb, 1e-6, None) / G_DAGGER))) * R)
        with np.errstate(invalid="ignore", divide="ignore"):
            n1 = (vo ** 2 - vb ** 2) / vo
        prof[gid] = dict(R=R, n1=n1, n2=2 * (vo - vb), vflat=sp.loc[gid, "Vflat"])
    return prof


def predict(prof, vc, r_kpc, key):
    v = [np.interp(r_kpc, p["R"], p[key]) for p in prof.values()
         if abs(p["vflat"] - vc) < 30 and p["R"].min() <= r_kpc <= p["R"].max()]
    return float(np.mean(v)) if v else np.nan


def manga_table():
    t = build()
    s = t[t.Vs_1 > 40]
    cr, pool = s[s.group == "CR"].copy(), s[s.group == "CO"]
    pairs = match(cr, pool, COLS)
    for c in (1.0, 1.5):
        for k in ("D_lo", "D_mid", "D_hi"):
            col = f"{k}_{c}"
            cr[f"s_{k}_{c}"] = [np.nanmedian(pool.loc[pairs[n], col]) - cr.loc[n, col]
                                if np.isfinite(pool.loc[pairs[n], col]).any() else np.nan
                                for n in cr.index]
    with ProcessPoolExecutor(8) as ex:
        cr["twist"] = list(ex.map(gas_twist, cr.index, cr.Re_arcsec))
    return cr


def compare(tab, prof, radius, label):
    out = {}
    for key in ("n1", "n2"):
        tab[f"pred_{key}"] = [predict(prof, g.Vc_star, radius * g.Re_kpc, key)
                              for _, g in tab.iterrows()]
        for k in ("D_mid", "D_lo", "D_hi"):
            d = pd.DataFrame({"twist": tab.twist, f"s_{k}_1.5": tab[f"s_{k}_{radius}"],
                              f"pred_{key}": tab[f"pred_{key}"]})
            d = d[np.isfinite(d[f"s_{k}_1.5"]) & np.isfinite(d.twist) & np.isfinite(d[f"pred_{key}"])]
            if len(d) >= 5:
                out[f"{key}_{k}"] = primary(d, float(d[f"pred_{key}"].mean()), k)
    r = out["n1_D_mid"]
    print(f"  {label} R = {radius} R_e: model prediction (psi = C0 r) {r['P']:.1f} km/s;"
          f" zero-twist slowdown a = {r['a']:+.1f} +/- {r['sd_a']:.1f};"
          f" f = {r['f']:.2f} (95% upper {r['f_upper95']:.2f}); log10 LR = {r['log10BF']:+.2f}"
          f" -> {r['label']} (N = {r['N']})")
    return out


def forecast(prof, cr):
    vc, re = float(cr.Vc_star.median()), float(cr.Re_kpc.median())
    rows = {}
    for m in (1.0, 1.5, 2.0, 3.0, 4.0):
        v = [np.interp(m * re, p["R"], p["n1"]) for p in prof.values()
             if abs(p["vflat"] - vc) < 30 and p["R"].min() <= m * re <= p["R"].max()]
        rows[m] = dict(pred=float(np.mean(v)) if v else np.nan, n_sparc=len(v))
    return dict(Vc_median=vc, Re_kpc_median=re, by_radius=rows)


def n_needed(P, sd_res, twist, n_sim=400):
    """Smallest N with P(BF > 10 for the true model) >= 0.8, both truths."""
    for N in (10, 20, 30, 40, 60, 80, 120):
        ok = []
        for truth, A in (("C", 0.0), ("V", P)):
            hits = 0
            for _ in range(n_sim):
                xs = RNG.choice(twist, N)
                ys = A + RNG.normal(0, sd_res, N)
                a, _, sd = intercept_and_sd(xs, ys, RNG, n_boot=200)
                lb = log10_bf(a, max(sd, 1e-6), P) * (1 if truth == "V" else -1)
                hits += lb > 1
            ok.append(hits / n_sim)
        if min(ok) >= 0.8:
            return N, ok
    return None, ok


def main():
    prof = model_profiles()
    res = {}
    print("POST HOC: counter-rotating gas vs the strong-gravitomagnetic GR model")
    cr = manga_table()
    res["MaNGA"] = {str(c): compare(cr, prof, c, "MaNGA") for c in (1.0, 1.5)}
    sami = pd.read_csv(ROOT / "results" / "c1_sample.csv", index_col=0)
    sami = sami[sami.role == "CR"]
    sami["pred_n1_15"] = [predict(prof, g.Vc_star, 1.5 * g.Re_kpc, "n1") for _, g in sami.iterrows()]
    sp = sami[np.isfinite(sami["s_D_mid_1.5"]) & np.isfinite(sami.twist)]
    res["SAMI_descriptive"] = dict(N=int(len(sp)), pred_mean=float(sp.pred_n1_15.mean()),
                                   s=sp["s_D_mid_1.5"].round(1).tolist(),
                                   twist=sp.twist.tolist())
    print(f"  SAMI (descriptive, N = {len(sp)}): model prediction {sp.pred_n1_15.mean():.1f} km/s;"
          f" s at 1.5 R_e = {sp['s_D_mid_1.5'].round(1).tolist()}")
    fc = forecast(prof, cr)
    res["forecast"] = fc
    print(f"  forecast for a typical counter-rotator (V_c = {fc['Vc_median']:.0f} km/s,"
          f" R_e = {fc['Re_kpc_median']:.1f} kpc):")
    for m, v in fc["by_radius"].items():
        print(f"    {m:.1f} R_e: predicted Delta V = {v['pred']:.0f} km/s ({v['n_sparc']} SPARC galaxies)")
    # residual scatter at 1.5 R_e from the MaNGA fit (robust), used for the outer-disc forecast
    d = cr[np.isfinite(cr["s_D_mid_1.5"]) & np.isfinite(cr.twist)]
    a, b, _ = intercept_and_sd(d.twist.values.astype(float), d["s_D_mid_1.5"].values, RNG)
    resid = d["s_D_mid_1.5"].values - (a + b * d.twist.values)
    sd_res = float(1.4826 * np.median(np.abs(resid - np.median(resid))))
    res["power_outer"] = {}
    for m in (2.0, 3.0):
        P = fc["by_radius"][m]["pred"]
        for infl in (1.0, 1.5):
            N, p = n_needed(P, sd_res * infl, d.twist.values.astype(float))
            res["power_outer"][f"{m}Re_noise_x{infl}"] = dict(P=P, sd=sd_res * infl, N_needed=N, power=p)
            print(f"    test at {m:.0f} R_e (P = {P:.0f}, scatter {sd_res * infl:.0f} km/s):"
                  f" N for 80% decisive power = {N}")
    cr.to_csv(ROOT / "results" / "agm_comparison_manga_cr.csv")
    (ROOT / "results" / "agm_comparison.json").write_text(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
