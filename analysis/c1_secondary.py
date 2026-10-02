#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""C1 secondary analyses 5, 7 and 8 (preregistration/OSF_C1_confirmatory.md).
Run after `analysis/c1_test.py run`. None of these decide the verdict.

  5  pooled MaNGA (exploratory) + C1 (confirmatory) primary statistic
  7  agreement of the pipeline's SAMI CR classification with the SAMI DR3
     catalogue PAs (PA_STELKIN, PA_GASKIN), and the primary statistic on
     the CR on which both agree
  8  the MaNGA primary statistic with and without the Galaxy Zoo merger cut

Usage: python analysis/c1_secondary.py {5,7,8}
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.c1_test import primary  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = ROOT / "results" / "c1_secondary.json"


def c1_primary_table():
    t = pd.read_csv(ROOT / "results" / "c1_sample.csv", index_col=0, dtype={0: str})
    cr = t[t.role == "CR"]
    return cr[np.isfinite(cr["s_D_mid_1.5"]) & np.isfinite(cr.twist)]


def sec5():
    m = pd.read_csv(ROOT / "results" / "g1_manga_discriminant_cr.csv", index_col=0)
    m = m[["s_15", "twist", "rz_15"]].dropna().rename(columns={"s_15": "s_D_mid_1.5"})
    c = c1_primary_table()[["s_D_mid_1.5", "twist", "rz_15"]]
    pool = pd.concat([m.assign(src="MaNGA"), c.assign(src="C1")])
    return {"pooled_MaNGA_C1": primary(pool, float(np.nanmean(pool.rz_15)), "D_mid"),
            "n_MaNGA": int(len(m)), "n_C1": int(len(c))}


def sec7():
    import pyvo
    c = c1_primary_table()
    s = c[c.survey == "SAMI"]
    tap = pyvo.dal.TAPService("https://datacentral.org.au/vo/tap")
    a = tap.run_sync("SELECT CATID, PA_STELKIN, PA_STELKIN_ERR FROM sami_dr3.samiDR3Stelkin",
                     maxrec=100000).to_table().to_pandas()
    b = tap.run_sync("SELECT CATID, PA_GASKIN, PA_GASKIN_ERR FROM sami_dr3.samiDR3gaskinPA",
                     maxrec=100000).to_table().to_pandas()
    p = a.merge(b, on="CATID").drop_duplicates("CATID")
    p["name"] = p.CATID.astype(str)
    p = p.set_index("name")
    d = (p.PA_GASKIN - p.PA_STELKIN).abs() % 360
    p["dpa_pub"] = np.minimum(d, 360 - d)
    j = s.join(p[["dpa_pub"]], how="left")
    agree = j[j.dpa_pub > 150]
    out = {"n_SAMI_CR": int(len(s)), "n_with_published_PA": int(j.dpa_pub.notna().sum()),
           "n_published_CR": int(len(agree))}
    keep = pd.concat([c[c.survey == "CALIFA"], agree.drop(columns="dpa_pub")])
    if len(keep) >= 10:
        out["primary_agreeing_CR"] = primary(keep, float(np.nanmean(keep.rz_15)), "D_mid")
    return out


def _manga_primary(with_merger_cut):
    from analysis.g1_manga_test import (gas_twist, load_catalogues, match,
                                        measure_one)
    from analysis.g1_manga_discriminant_exploratory import pred_rz
    from analysis.s_pred_rar import rar_profiles
    d = load_catalogues()
    d = d[(d.e_pa_star <= 20) & (d.e_pa_gas <= 20) & d.dpa.notna()]
    if with_merger_cut:
        d = d[~d.merger]
    d["group"] = np.where(d.dpa > 150, "CR", np.where(d.dpa < 30, "CO", ""))
    d = d[d.group != ""].copy()
    d["Re_kpc"] = d.Re_arcsec * d.kpc_per_arcsec
    d["log_sig1re"] = np.log10(d.sig1re.clip(lower=1))
    jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in d.iterrows()]
    rec = []
    with ProcessPoolExecutor(8) as ex:
        for pifu, rows, _ in ex.map(measure_one, jobs, chunksize=4):
            if rows is None:
                continue
            r = {q["r_re"]: q for q in rows}
            rec.append({"plateifu": pifu, "Vs_1": r[1.0]["Vs"], "Vc_star": r[1.0]["Vc_star"],
                        "D15": r[1.5]["D_mid"]})
    t = d.join(pd.DataFrame(rec).set_index("plateifu"), how="inner")
    cols = ["logMstar", "ttype", "log_sig1re", "inc"]
    t = t[np.isfinite(t.Vs_1) & (t.Vs_1 > 40)].dropna(subset=cols)
    cr, pool = t[t.group == "CR"], t[t.group == "CO"]
    pairs = match(cr, pool, cols)
    cr = cr.copy()
    cr["s_D_mid_1.5"] = [np.nanmedian(pool.loc[pairs[n], "D15"]) - cr.loc[n, "D15"]
                         if np.isfinite(pool.loc[pairs[n], "D15"]).any() else np.nan
                         for n in cr.index]
    with ProcessPoolExecutor(8) as ex:
        cr["twist"] = list(ex.map(gas_twist, cr.index, cr.Re_arcsec))
    prof = rar_profiles()
    cr["rz_15"] = [pred_rz(prof, g.Vc_star, 1.5 * g.Re_kpc) for _, g in cr.iterrows()]
    prim = cr[np.isfinite(cr["s_D_mid_1.5"]) & np.isfinite(cr.twist)]
    return primary(prim, float(np.nanmean(prim.rz_15)), "D_mid")


def sec8():
    return {"MaNGA_with_merger_cut": _manga_primary(True),
            "MaNGA_without_merger_cut": _manga_primary(False)}


def main():
    which = sys.argv[1]
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    res[which] = {"5": sec5, "7": sec7, "8": sec8}[which]()
    print(json.dumps(res[which], indent=1, default=float))
    OUT.write_text(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
