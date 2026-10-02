#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Confirmatory test C1 (SAMI DR3 + CALIFA): preregistration/OSF_C1_confirmatory.md.

Stages (run in this order):
  select   kinematic PAs for the parent sample -> data/processed/c1_pa.csv
           (PAs, errors, spaxel counts and Delta PA only)
  dry      groups, measurements, cuts, matching and twist; prints sample
           counts only (no D, s, a, b, BF or f)
  run      the registered analysis, once -> results/c1.json, results/c1_sample.csv

Usage: python analysis/c1_test.py {select,dry,run} [--workers W] [--limit N]
"""

import argparse
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import califa, sami, manga, manga_measure as mm  # noqa: E402
from analysis.g1_manga_select import receding_pa  # noqa: E402
from analysis.g1_manga_test import boot_S, match  # noqa: E402
from analysis.g1_manga_discriminant_exploratory import pred_rz  # noqa: E402
from analysis.s_pred_rar import rar_profiles  # noqa: E402
from analysis.c1_power import theil_sen  # noqa: E402

warnings.filterwarnings("ignore")
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SEED = 20261002
N_BOOT = 10_000
V_CUT = 40.0
OVERLAP_ARCSEC = 5.0
PA_MAX_ERR = 20.0
MIN_SPAX = 50                # valid spaxels after 2x decimation, as in MaNGA
DECIMATE = 2
COVS = ["logMstar", "morph", "log_sig1re", "inc"]
PA_OUT = ROOT / "data" / "processed" / "c1_pa.csv"
MEAS_OUT = ROOT / "data" / "processed" / "c1_meas.csv"
MOD = {"SAMI": sami, "CALIFA": califa}


# ---------------------------------------------------------------- sample
def parent_sample():
    s, c = sami.catalogue(), califa.catalogue()
    keep = ["ra", "dec", "z", "Re_arcsec", "inc", "logMstar", "morph",
            "kpc_per_arcsec", "survey", "quality_ok"]
    d = pd.concat([c[keep], s[keep]])
    log = {"parent_CALIFA": int(len(c)), "parent_SAMI": int(len(s))}
    d = d[d.quality_ok & d.inc.between(30, 80) & (d.Re_arcsec > 0)
          & np.isfinite(d.kpc_per_arcsec)]
    # drop anything within 5" of a MaNGA DR17 target
    drp = fits.getdata(manga.MANGA / "drpall-v3_1_1.fits", "MANGA")
    mc = SkyCoord(drp["objra"].astype(float) * u.deg, drp["objdec"].astype(float) * u.deg)
    pc = SkyCoord(d.ra.values * u.deg, d.dec.values * u.deg)
    _, sep, _ = pc.match_to_catalog_sky(mc)
    in_manga = sep.arcsec < OVERLAP_ARCSEC
    log["dropped_in_MaNGA"] = int(in_manga.sum())
    d = d[~in_manga]
    # SAMI galaxies also in CALIFA: keep the CALIFA entry
    cs, ss = d[d.survey == "CALIFA"], d[d.survey == "SAMI"]
    if len(cs) and len(ss):
        _, sep, _ = SkyCoord(ss.ra.values * u.deg, ss.dec.values * u.deg).match_to_catalog_sky(
            SkyCoord(cs.ra.values * u.deg, cs.dec.values * u.deg))
        dup = ss.index[sep.arcsec < OVERLAP_ARCSEC]
        log["dropped_SAMI_in_CALIFA"] = int(len(dup))
        d = d.drop(dup)
    log["parent_after_cuts"] = int(len(d))
    return d, log


def arrays(name, row):
    if row.survey == "CALIFA":
        return califa.galaxy_arrays(name, row.ra, row.dec, row.z)
    return sami.galaxy_arrays(name)


# ---------------------------------------------------------------- stage 1
def pa_one(args):
    name, row = args
    maps, x, y, s_ok, g_ok = arrays(name, row)
    dec = np.zeros_like(s_ok)
    dec[::DECIMATE, ::DECIMATE] = True
    s, g = s_ok & dec, g_ok & dec
    out = {"name": name, "survey": row.survey, "n_star": int(s.sum()), "n_gas": int(g.sum())}
    if out["n_star"] >= MIN_SPAX:
        out["pa_star"], out["e_pa_star"] = receding_pa(x[s], y[s], maps["vs"][s], maps["evs"][s])
    if out["n_gas"] >= MIN_SPAX:
        out["pa_gas"], out["e_pa_gas"] = receding_pa(x[g], y[g], maps["vg"][g], maps["evg"][g])
    if "pa_star" in out and "pa_gas" in out:
        dd = abs(out["pa_gas"] - out["pa_star"]) % 360.0
        out["dpa"] = min(dd, 360.0 - dd)
    return out


def stage_select(args):
    d, log = parent_sample()
    sami_missing = [n for n, r in d.iterrows() if r.survey == "SAMI" and not sami.available(n)]
    log["SAMI_maps_missing"] = len(sami_missing)
    d = d.drop(sami_missing)
    print(json.dumps(log, indent=1), flush=True)
    done = set(pd.read_csv(PA_OUT)["name"].astype(str)) if PA_OUT.exists() else set()
    todo = [(n, r) for n, r in d.iterrows() if n not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"PA fits: done {len(done)}, to do {len(todo)}", flush=True)
    rows, fails = [], 0
    with ProcessPoolExecutor(args.workers) as ex:
        futs = {ex.submit(pa_one, t): t[0] for t in todo}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                rows.append(f.result())
            except Exception as e:  # noqa: BLE001  (download/map error: retried next run)
                fails += 1
                print(f"  fail {futs[f]}: {type(e).__name__}: {e}", flush=True)
            if rows and (i % 25 == 0 or i == len(todo)):
                pd.DataFrame(rows).to_csv(PA_OUT, mode="a", index=False, header=not PA_OUT.exists())
                rows = []
                print(f"  {i}/{len(todo)} ({fails} failed)", flush=True)


# ---------------------------------------------------------------- stage 2
def twist(maps, x, y, re):
    r = np.hypot(x, y)
    g = np.isfinite(maps["vg"])
    out = []
    for sel in (g & (r < re), g & (r >= re)):
        if sel.sum() < 30:
            return np.nan
        out.append(receding_pa(x[sel], y[sel], maps["vg"][sel], maps["evg"][sel])[0])
    dd = abs(out[0] - out[1]) % 360
    return min(dd, 360 - dd)


def measure_one(args):
    name, row = args
    try:
        maps, x, y, _, g_ok = arrays(name, row)
        rows = mm.measure(maps, x, y, row.pa_star, row.pa_gas, row.inc, row.Re_arcsec)
        R, _ = mm.disk_coords(x, y, row.pa_star, row.inc)
        inside = (R < row.Re_arcsec) & np.isfinite(maps["ss"])
        out = {"name": name,
               "sig1re": float(np.median(maps["ss"][inside])) if inside.sum() >= 5 else np.nan,
               "med_eg": float(np.median(maps["evg"][g_ok])) if g_ok.any() else np.nan,
               "twist": twist(maps, x, y, row.Re_arcsec)}
        r = {q["r_re"]: q for q in rows}
        for c in (1.0, 1.5):
            for k in ("D_lo", "D_mid", "D_hi"):
                out[f"{k}_{c}"] = r[c][k]
        out["Vs_1"], out["Vc_star"] = r[1.0]["Vs"], r[1.0]["Vc_star"]
        out["sg_15"] = r[1.5]["sg"]
        out["gas_resid"] = np.nanmean([q["gas_resid"] for q in rows])
        return out
    except Exception as e:  # noqa: BLE001
        return {"name": name, "error": repr(e)}


def build(args):
    d, log = parent_sample()
    pa = pd.read_csv(PA_OUT, dtype={"name": str}).drop_duplicates("name", keep="last")
    pa = pa.set_index("name").drop(columns="survey")
    d = d.join(pa, how="inner")
    ok = (d.e_pa_star <= PA_MAX_ERR) & (d.e_pa_gas <= PA_MAX_ERR) & d.dpa.notna()
    d = d[ok].copy()
    d["group"] = np.where(d.dpa > 150, "CR", np.where(d.dpa < 30, "CO", ""))
    d = d[d.group != ""]
    log["CR_candidates"] = int((d.group == "CR").sum())
    log["CO_candidates"] = int((d.group == "CO").sum())
    if MEAS_OUT.exists():
        meas = pd.read_csv(MEAS_OUT, dtype={"name": str}).set_index("name")
    else:
        with ProcessPoolExecutor(args.workers) as ex:
            meas = pd.DataFrame(list(ex.map(measure_one, d.iterrows(), chunksize=2))).set_index("name")
        meas.to_csv(MEAS_OUT)
    if "error" in meas:
        log["measure_failed"] = int(meas["error"].notna().sum())
        meas = meas[meas["error"].isna()]
    t = d.join(meas, how="inner")
    t["Re_kpc"] = t.Re_arcsec * t.kpc_per_arcsec
    t["log_sig1re"] = np.log10(t.sig1re.clip(lower=1))
    t = t[np.isfinite(t.Vs_1) & (t.Vs_1 > V_CUT)].dropna(subset=COVS)
    log["CR_after_cuts"] = int((t.group == "CR").sum())
    log["CO_pool_after_cuts"] = int((t.group == "CO").sum())
    # match within survey
    cr_all, co_ids, ctrl = [], [], {}
    for sv in ("SAMI", "CALIFA"):
        ts = t[t.survey == sv]
        cr, pool = ts[ts.group == "CR"], ts[ts.group == "CO"]
        log[f"CR_{sv}"] = int(len(cr))
        if len(cr) == 0:
            continue
        n_per = controls_per_cr(len(cr), len(pool))
        log[f"controls_per_CR_{sv}"] = n_per
        if n_per == 0:
            continue
        pairs = match(cr, pool, COVS, n_per=n_per)
        ctrl.update(pairs)
        cr_all.append(cr)
        co_ids += [p for v in pairs.values() for p in v]
    cr = pd.concat(cr_all) if cr_all else t.iloc[:0]
    co = t.loc[co_ids]
    for c in (1.0, 1.5):
        for k in ("D_lo", "D_mid", "D_hi"):
            cr[f"s_{k}_{c}"] = [np.nanmedian(co.loc[ctrl[n], f"{k}_{c}"]) - cr.loc[n, f"{k}_{c}"]
                                if np.isfinite(co.loc[ctrl[n], f"{k}_{c}"]).any() else np.nan
                                for n in cr.index]
    prim = cr[np.isfinite(cr["s_D_mid_1.5"]) & np.isfinite(cr.twist)]
    log["CR_matched"] = int(len(cr))
    log["CO_matched"] = int(len(co))
    log["CR_primary_1.5Re"] = int(len(prim))
    for sv in ("SAMI", "CALIFA"):
        log[f"CR_primary_{sv}"] = int((prim.survey == sv).sum())
    return t, cr, co, prim, log


def controls_per_cr(n_cr, n_pool):
    """3 controls per CR if the pool allows, else 2, else 1, else 0 (excluded)."""
    for n in (3, 2, 1):
        if n_pool >= n * n_cr:
            return n
    return 0


def balance(cr, co):
    """Standardised mean differences (CR - CO) / pooled SD; reported only."""
    out = {}
    for c in COVS + ["Vc_star", "sg_15"]:
        a, b = cr[c].dropna(), co[c].dropna()
        sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        out[c] = dict(mean_CR=float(a.mean()), mean_CO=float(b.mean()),
                      smd=float((a.mean() - b.mean()) / sd) if sd > 0 else np.nan)
    return out


# ---------------------------------------------------------------- stage 3
def intercept_boot(x, y, seed=SEED, n=N_BOOT):
    a, b = theil_sen(x, y)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), (n, len(x)))
    ab, bb = [], []
    for chunk in np.array_split(idx, 20):
        aa, bbb = theil_sen(x[chunk], y[chunk])
        ab.append(aa)
        bb.append(bbb)
    return float(a[0]), float(b[0]), np.concatenate(ab), np.concatenate(bb)


def label(bf):
    if bf >= 10:
        return 2, "strong evidence for H_V"
    if bf >= 3:
        return 1, "moderate evidence for H_V"
    if bf <= 0.1:
        return -2, "strong evidence for H_C"
    if bf <= 1 / 3:
        return -1, "moderate evidence for H_C"
    return 0, "inconclusive"


def primary(prim, P, k):
    x, y = prim.twist.values.astype(float), prim[f"s_{k}_1.5"].values
    a, b, ab, bb = intercept_boot(x, y)
    sd = float(np.nanstd(ab))
    log10bf = (a ** 2 - (a - P) ** 2) / (2 * sd ** 2) / np.log(10)
    f_lo, f_hi = np.nanpercentile(ab / P, [2.5, 97.5])
    return dict(a=a, sd_a=sd, b=b, b_lo95_onesided=float(np.nanpercentile(bb, 5)),
                P=P, log10BF=float(log10bf), BF=float(10 ** log10bf),
                f=a / P, f_ci95=[float(f_lo), float(f_hi)], f_upper95=float(np.nanpercentile(ab / P, 95)),
                label=label(10 ** log10bf)[1], code=label(10 ** log10bf)[0], N=int(len(prim)))


def stage_run(args):
    t, cr, co, prim, log = build(args)
    print(json.dumps(log, indent=1), flush=True)
    if args.stage == "dry":
        print("DRY RUN: mechanics OK; no D, s, a, b, BF or f computed or printed.")
        return
    prof = rar_profiles()
    cr["rz_15"] = [pred_rz(prof, g.Vc_star, 1.5 * g.Re_kpc) for _, g in cr.iterrows()]
    cr["rz_1"] = [pred_rz(prof, g.Vc_star, 1.0 * g.Re_kpc) for _, g in cr.iterrows()]
    prim = cr.loc[prim.index]
    res = {"log": log}
    if len(prim) < 10:
        res["verdict"] = f"UNDERPOWERED (N = {len(prim)} < 10): descriptive only"
        P = float(np.nanmean(prim.rz_15)) if len(prim) else np.nan
        res["descriptive"] = dict(P=P, s_median=float(np.nanmedian(prim["s_D_mid_1.5"])) if len(prim) else None)
    else:
        P = float(np.nanmean(prim.rz_15))
        res["primary"] = {k: primary(prim, P, k) for k in ("D_mid", "D_lo", "D_hi")}
        codes = [res["primary"][k]["code"] for k in ("D_mid", "D_lo", "D_hi")]
        if len(set(np.sign(codes))) > 1:
            head = 0
        else:
            head = int(np.sign(codes[0]) * min(abs(c) for c in codes))
        names = {2: "strong evidence for H_V", 1: "moderate evidence for H_V", 0: "inconclusive",
                 -1: "moderate evidence for H_C", -2: "strong evidence for H_C"}
        res["verdict"] = names[head]
        # secondary analyses (preregistered; not decisive)
        sec = {}
        for c in (1.0, 1.5):  # 1: S as in MaNGA
            for k in ("D_lo", "D_mid", "D_hi"):
                a_ = cr[f"{k}_{c}"].dropna().values
                b_ = co[f"{k}_{c}"].dropna().values
                if len(a_) > 2 and len(b_) > 2:
                    S, sS = boot_S(a_, b_)
                    sec[f"S_{k}_{c}"] = dict(S=S, sS=sS, n_cr=len(a_), n_co=len(b_))
        for sv in ("SAMI", "CALIFA"):  # 4: per survey
            ps = prim[prim.survey == sv]
            if len(ps) >= 5:
                sec[f"primary_{sv}"] = primary(ps, float(np.nanmean(ps.rz_15)), "D_mid")
        from scipy import stats  # 6: Spearman with predicted R*zeta
        m = np.isfinite(prim.rz_15)
        rho = stats.spearmanr(prim.rz_15[m], prim["s_D_mid_1.5"][m])
        sec["spearman_s_rz"] = dict(rho=float(rho.statistic), p=float(rho.pvalue), n=int(m.sum()))
        sec["balance"] = balance(cr, co)  # 9: CR-control balance
        res["secondary"] = sec
        print(json.dumps(res, indent=1, default=float))
    print(f"\nREGISTERED VERDICT (C1): {res['verdict']}")
    pd.concat([cr.assign(role="CR"), co.assign(role="CO")]).to_csv(ROOT / "results" / "c1_sample.csv")
    (ROOT / "results" / "c1.json").write_text(json.dumps(res, indent=2, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("select", "dry", "run"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    if args.stage == "select":
        stage_select(args)
    else:
        stage_run(args)


if __name__ == "__main__":
    main()
