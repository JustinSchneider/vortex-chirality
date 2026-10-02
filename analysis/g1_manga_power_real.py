#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Realistic power check for MaNGA G1, using the real sample's nuisance
properties but never the test statistic.

Allowed inputs per galaxy: stellar V and sigma, gas sigma, drift factors k,
ring-fit velocity errors, gas-disturbance index, R_e, matching variables.
Forbidden: gas rotation speed and D. measure() returns those too; they are
dropped immediately and never printed or saved.

For each candidate speed cut, the real counter-rotators (CR) and their
registered matches (CO) are used. Synthetic D values are then simulated:
  D_i = s_i R*zeta_i (CR only, signal runs) + gas scatter + per-galaxy k error
        + measurement error + shared k offset (unknown truth, applied as in
          the registered k-1/k/k+1 rule)
and the registered verdict rules are applied.
"""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import manga_measure as mm, models  # noqa: E402
from analysis.g1_manga_test import galaxy_arrays, load_catalogues, match  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NUIS = ROOT / "data" / "processed" / "manga_nuisance.csv"
ALLOWED = ("Vs", "eVs", "ss", "ks", "eVg", "sg", "kg", "gas_resid", "Vc_star")
RNG = np.random.default_rng(20261002)
N_SIM = 300
N_BOOT = 300
GAS_FRAC = 0.15
SIG_K = 0.5


def nuisance_one(args):
    pifu, pa_s, pa_g, inc, re = args
    try:
        maps, x, y, _, g_ok = galaxy_arrays(pifu)
        rows = mm.measure(maps, x, y, pa_s, pa_g, inc, re)
        r1 = next(r for r in rows if r["r_re"] == 1.0)
        out = {k: r1[k] for k in ALLOWED}          # drops Vg and all D values
        out["gas_resid"] = np.nanmean([r["gas_resid"] for r in rows])
        out["med_eg"] = float(np.median(maps["evg"][g_ok])) if g_ok.any() else np.nan
        out["plateifu"] = pifu
        return out
    except Exception:  # noqa: BLE001
        return None


def build_nuisance(d):
    jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in d.iterrows()]
    with ProcessPoolExecutor(12) as ex:
        rows = [r for r in ex.map(nuisance_one, jobs, chunksize=4) if r]
    t = pd.DataFrame(rows).set_index("plateifu")
    assert not {"Vg", "D_lo", "D_mid", "D_hi"} & set(t.columns)
    t.to_csv(NUIS)
    return t


def sparc_pool():
    c = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv", index_col=0)
    return c[c["clean"] & (c["Vflat"] > 0)]


def rzeta_draw(sp, Vc, R_kpc):
    ref = sp[(sp["Vflat"] - Vc).abs() < 30]
    if len(ref) == 0:
        return np.array([0.0])
    return R_kpc * models.rt_flow_vorticity(R_kpc, ref["omega"].values, ref["Rt"].values)


def boot_E(a, b):
    E = a.var(ddof=1) - b.var(ddof=1)
    bs = [RNG.choice(a, len(a)).var(ddof=1) - RNG.choice(b, len(b)).var(ddof=1)
          for _ in range(N_BOOT)]
    return E, np.std(bs, ddof=1)


def simulate(cr, co, f_dist, sp, signal):
    """One synthetic experiment; returns dict of registered verdict flags."""
    def nuis(g, f):
        Vc = g["Vc_star"].values
        gain = (g["sg"].fillna(0).values ** 2 - g["ss"].values ** 2) / (2 * Vc)
        noise = (f * RNG.normal(0, GAS_FRAC * Vc)
                 + RNG.normal(0, SIG_K, len(g)) * gain
                 + RNG.normal(0, np.sqrt(g["eVs"].fillna(5) ** 2
                                         + g["eVg"].fillna(5) ** 2)))
        return noise, gain

    n_cr, g_cr = nuis(cr, f_dist)
    n_co, g_co = nuis(co, 1.0)
    sig = np.zeros(len(cr))
    if signal:
        for i, (_, g) in enumerate(cr.iterrows()):
            sig[i] = RNG.choice([-1, 1]) * RNG.choice(rzeta_draw(sp, g["Vc_star"], g["Re_kpc"]))
    delta = RNG.uniform(-1, 1)                      # unknown true k offset
    res = []
    for j in (-1, 0, 1):
        a = sig + n_cr + (j - delta) * g_cr
        b = n_co + (j - delta) * g_co
        res.append(boot_E(a, b))
    return res


def verdicts(res, E_pred):
    fals = all(E + 2 * s < 0.25 * E_pred for E, s in res)
    supp = all(E > 3 * s and abs(E - E_pred) < 2 * s for E, s in res)
    sig = all(E > 3 * s for E, s in res)
    return fals, supp, sig


def main():
    d = load_catalogues()
    d = d[(d["e_pa_star"] <= 20) & (d["e_pa_gas"] <= 20) & d["dpa"].notna() & ~d["merger"]]
    d["group"] = np.where(d["dpa"] > 150, "CR", np.where(d["dpa"] < 30, "CO", ""))
    d = d[d["group"] != ""].copy()
    d["Re_kpc"] = d["Re_arcsec"] * d["kpc_per_arcsec"]
    d["log_sig1re"] = np.log10(d["sig1re"].clip(lower=1))
    t = pd.read_csv(NUIS, index_col=0) if NUIS.exists() else build_nuisance(d)
    t = d.join(t, how="inner")
    t = t[np.isfinite(t["Vs"]) & np.isfinite(t["ss"]) & (t["Vc_star"] > 0)]
    t = t.dropna(subset=["logMstar", "ttype", "log_sig1re", "inc"])
    t["resid_ratio"] = t["gas_resid"] / t["med_eg"]
    sp = sparc_pool()

    out = {}
    print("=" * 78)
    print("Realistic power check (real nuisance properties; no gas speeds or D used)")
    print("=" * 78)
    for vcut in (100, 80, 60, 40):
        s = t[t["Vs"] > vcut]
        cr, pool = s[s.group == "CR"], s[s.group == "CO"]
        if len(cr) < 5 or len(pool) < 3 * len(cr):
            continue
        pairs = match(cr, pool, ["logMstar", "ttype", "log_sig1re", "inc"])
        co = pool.loc[[p for v in pairs.values() for p in v]]
        f_dist = float(np.nanmedian(cr["resid_ratio"]) / np.nanmedian(co["resid_ratio"]))
        E_pred = float(np.mean([np.mean(rzeta_draw(sp, g["Vc_star"], g["Re_kpc"]) ** 2)
                                for _, g in cr.iterrows()]))
        gain = (cr["sg"].fillna(0) ** 2 - cr["ss"] ** 2).abs() / (2 * cr["Vc_star"])
        print(f"\nV_* > {vcut} km/s: N_CR = {len(cr)}, N_CO = {len(co)}")
        print(f"  CR median V_* {cr['Vs'].median():.0f}, sigma_* {cr['ss'].median():.0f} km/s;"
              f" drift sensitivity |dD/dk| median {gain.median():.0f} km/s per unit k")
        print(f"  disturbance ratio CR/CO (gas residual) f = {f_dist:.2f};"
              f"  E_pred = {E_pred:.0f} (RMS R*zeta {np.sqrt(E_pred):.0f} km/s)")
        row = {"N_CR": len(cr), "f_dist": f_dist, "E_pred": E_pred}
        for f_lab, f in (("measured f", f_dist), ("f x 1.3", f_dist * 1.3)):
            cnt = np.zeros((2, 3))
            for _ in range(N_SIM):
                for si, signal in enumerate((True, False)):
                    cnt[si] += verdicts(simulate(cr, co, f, sp, signal), E_pred)
            p = cnt / N_SIM
            row[f_lab] = dict(supp_if_signal=p[0, 1], fals_if_signal=p[0, 0],
                              fals_if_null=p[1, 0], supp_if_null=p[1, 1],
                              sig_if_null=p[1, 2])
            print(f"  [{f_lab:10s}] signal present: SUPPORTED {p[0,1]:.2f}, falsified {p[0,0]:.2f}"
                  f" | signal absent: FALSIFIED {p[1,0]:.2f}, false-support {p[1,1]:.2f},"
                  f" E>3σ {p[1,2]:.2f}")
        out[str(vcut)] = row
    (ROOT / "results" / "g1_manga_power_real.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
