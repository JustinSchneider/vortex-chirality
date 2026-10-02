#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""EXPLORATORY: is the counter-rotator slowdown vortex-like or
disequilibrium-like? (Informs the confirmatory pre-registration; not a test.)

Per counter-rotator i and radius R in {1, 1.5} R_e (central k):
  slowdown s_i = median(D_CO) - D_CR,i
Vortex-like (velocity-dependent field aligned with the stars):
  s_i tracks the galaxy's predicted R*zeta_i (RAR normalisation, evaluated
  at its own radius), and s(1.5)/s(1) follows R*zeta(1.5)/R*zeta(1).
Disequilibrium-like (accreted gas not yet settled):
  s_i tracks gas disturbance: residual ratio, kinematic twist, |V_R|, and
  pressure support sigma_gas / V_c.
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
from analysis.s_pred_rar import rar_profiles  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
RNG = np.random.default_rng(20261005)


def extra(args):
    pifu, pa_s, pa_g, inc, re = args
    maps, x, y, _, _ = galaxy_arrays(pifu)
    r = {row["r_re"]: row for row in mm.measure(maps, x, y, pa_s, pa_g, inc, re)}
    return pifu, {"sg_1": r[1.0]["sg"], "VR_1": r[1.0]["VRg"],
                  "sg_15": r[1.5]["sg"], "VR_15": r[1.5]["VRg"]}


def pred_rz(prof, Vc, R_kpc):
    rz = [np.interp(R_kpc, R, z) for R, z, vf in prof.values()
          if abs(vf - Vc) < 30 and R.min() <= R_kpc <= R.max()]
    return float(np.mean(rz)) if rz else np.nan


def spear(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    r = stats.spearmanr(x[m], y[m])
    return float(r.statistic), float(r.pvalue), int(m.sum())


def main():
    t = pd.read_csv(ROOT / "results" / "g1_manga_sample.csv", index_col=0)
    cr, co = t[t.group == "CR"].copy(), t[t.group == "CO"]
    jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in cr.iterrows()]
    with ProcessPoolExecutor(10) as ex:
        ext = dict(ex.map(extra, jobs))
    for k in ("sg_1", "VR_1", "sg_15", "VR_15"):
        cr[k] = [ext[p][k] for p in cr.index]
    prof = rar_profiles()
    out = {}
    print("EXPLORATORY discriminant analysis (MaNGA, V_* > 40, central k)")
    for rad, tag in ((1.0, "1"), (1.5, "15")):
        s = co[f"D_mid_{rad}"].median() - cr[f"D_mid_{rad}"]
        cr[f"s_{tag}"] = s
        cr[f"rz_{tag}"] = [pred_rz(prof, g.Vc_star, rad * g.Re_kpc) for _, g in cr.iterrows()]
        res = {}
        print(f"\n  R = {rad} R_e: median slowdown {np.nanmedian(s):+.1f} km/s,"
              f" median predicted R*zeta {np.nanmedian(cr[f'rz_{tag}']):.1f} km/s")
        tests = {
            "VORTEX: predicted R*zeta": cr[f"rz_{tag}"],
            "DISEQ: sigma_gas / V_c": cr[f"sg_{tag}"] / cr["Vc_star"],
            "DISEQ: gas residual ratio": cr["resid_ratio"],
            "DISEQ: gas PA twist": cr["twist"],
            "DISEQ: |V_R| / V_c": cr[f"VR_{tag}"] / cr["Vc_star"],
            "(context) V_c,*": cr["Vc_star"],
        }
        for name, x in tests.items():
            rho, p, n = spear(x.values, s.values)
            res[name] = dict(rho=rho, p=p, n=n)
            print(f"    rho(slowdown, {name:28s}) = {rho:+.2f}  (p = {p:.3f}, N = {n})")
        # vortex: partial correlation with predicted R*zeta at fixed V_c
        m = np.isfinite(s) & np.isfinite(cr[f"rz_{tag}"])
        rk = lambda a: stats.rankdata(a)
        zz = rk(cr["Vc_star"][m])
        resid = lambda a: a - np.polyval(np.polyfit(zz, a, 1), zz)
        pr = stats.pearsonr(resid(rk(s[m])), resid(rk(cr[f"rz_{tag}"][m])))
        res["partial_rz_given_Vc"] = dict(rho=float(pr.statistic), p=float(pr.pvalue))
        print(f"    partial rho(slowdown, R*zeta | V_c)          = {pr.statistic:+.2f}  (p = {pr.pvalue:.3f})")
        out[str(rad)] = res
    # radial growth: observed vs vortex-predicted
    both = cr[np.isfinite(cr.s_1) & np.isfinite(cr.s_15) & (cr.rz_1 > 0)]
    obs_ratio = both.s_15.median() / both.s_1.median() if both.s_1.median() != 0 else np.nan
    pred_ratio = (both.rz_15 / both.rz_1).median()
    boots = []
    for _ in range(5000):
        b = both.sample(len(both), replace=True, random_state=int(RNG.integers(1 << 31)))
        boots.append(b.s_15.median() - b.s_1.median())
    print(f"\n  radial growth: median slowdown 1 R_e {both.s_1.median():+.1f}, 1.5 R_e {both.s_15.median():+.1f}"
          f" (difference {both.s_15.median() - both.s_1.median():+.1f} ± {np.std(boots):.1f});"
          f" vortex-predicted R*zeta ratio 1.5/1 = {pred_ratio:.2f}")
    out["radial"] = dict(obs_s1=float(both.s_1.median()), obs_s15=float(both.s_15.median()),
                         diff_sd=float(np.std(boots)), pred_ratio=float(pred_ratio), n=len(both))
    cr.to_csv(ROOT / "results" / "g1_manga_discriminant_cr.csv")
    (ROOT / "results" / "g1_manga_discriminant_exploratory.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
