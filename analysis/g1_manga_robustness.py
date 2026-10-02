#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Robustness of the MaNGA chirality result to analysis choices (post hoc).

Grid: stellar speed cut V_*(1 R_e) in {40, 60, 80, 100} km/s  x  radius
{1.0, 1.5} R_e  x  drift factor {k-1, k, k+1}  x  amplitude normalisation
{flow-profile template (pre-specified), RAR}. Matching and statistics reuse
analysis/g1_manga_test.py unchanged. Measurements for all candidates are
cached in data/processed/manga_meas_all.csv.
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.g1_manga_test import (boot_S, load_catalogues, match,  # noqa: E402
                                    measure_one, s_pred)
from analysis.s_pred_rar import rar_profiles  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
CACHE = ROOT / "data" / "processed" / "manga_meas_all.csv"
COLS = ["logMstar", "ttype", "log_sig1re", "inc"]


def build():
    d = load_catalogues()
    d = d[(d.e_pa_star <= 20) & (d.e_pa_gas <= 20) & d.dpa.notna() & ~d.merger]
    d["group"] = np.where(d.dpa > 150, "CR", np.where(d.dpa < 30, "CO", ""))
    d = d[d.group != ""].copy()
    d["Re_kpc"] = d.Re_arcsec * d.kpc_per_arcsec
    d["log_sig1re"] = np.log10(d.sig1re.clip(lower=1))
    if CACHE.exists():
        meas = pd.read_csv(CACHE, index_col=0)
    else:
        jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in d.iterrows()]
        rec = []
        with ProcessPoolExecutor(10) as ex:
            for pifu, rows, _ in ex.map(measure_one, jobs, chunksize=4):
                if rows is None:
                    continue
                r = {row["r_re"]: row for row in rows}
                out = {"plateifu": pifu, "Vs_1": r[1.0]["Vs"], "Vc_star": r[1.0]["Vc_star"]}
                for c in (1.0, 1.5):
                    for k in ("D_lo", "D_mid", "D_hi"):
                        out[f"{k}_{c}"] = r[c][k]
                rec.append(out)
        meas = pd.DataFrame(rec).set_index("plateifu")
        meas.to_csv(CACHE)
    t = d.join(meas, how="inner").dropna(subset=COLS)
    return t[np.isfinite(t["D_mid_1.0"])]


def s_pred_rar_for(cr, prof):
    vals = []
    for _, g in cr.iterrows():
        rz = [np.interp(g.Re_kpc, R, z) for R, z, vf in prof.values()
              if abs(vf - g.Vc_star) < 30 and R.min() <= g.Re_kpc <= R.max()]
        if rz:
            vals.append(np.mean(rz))
    return float(np.mean(vals))


def main():
    t = build()
    prof = rar_profiles()
    rows = []
    for vcut in (40, 60, 80, 100):
        s = t[t.Vs_1 > vcut]
        cr, pool = s[s.group == "CR"], s[s.group == "CO"]
        if len(cr) < 10 or len(pool) < 3 * len(cr):
            continue
        pairs = match(cr, pool, COLS)
        co = pool.loc[[p for v in pairs.values() for p in v]]
        sp_t, _ = s_pred(cr)
        sp_r = s_pred_rar_for(cr, prof)
        for c in (1.0, 1.5):
            for k, lab in (("D_lo", "k-1"), ("D_mid", "k"), ("D_hi", "k+1")):
                a = cr[f"{k}_{c}"].dropna().values
                b = co[f"{k}_{c}"].dropna().values
                S, sS = boot_S(a, b)
                rows.append(dict(V_cut=vcut, R_Re=c, k=lab, N_CR=len(a), N_CO=len(b),
                                 S=S, sS=sS, S_pred_template=sp_t, S_pred_RAR=sp_r,
                                 excl_template_sigma=(sp_t - S) / sS,
                                 excl_RAR_sigma=(sp_r - S) / sS,
                                 upper2s_frac_template=(S + 2 * sS) / sp_t,
                                 upper2s_frac_RAR=(S + 2 * sS) / sp_r))
    r = pd.DataFrame(rows)
    r.to_csv(ROOT / "results" / "g1_manga_robustness.csv", index=False)
    pd.set_option("display.width", 200)
    print(r.round(2).to_string(index=False))
    # Note: S_pred is evaluated at R_e for both radii (as pre-specified); at
    # 1.5 R_e the comparison is indicative only.


if __name__ == "__main__":
    main()
