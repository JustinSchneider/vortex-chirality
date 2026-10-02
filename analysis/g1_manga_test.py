#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""MaNGA G1, stage 2: the registered chirality test (PREREG amendments c-f).

Steps: groups from Delta PA -> merger exclusion -> ring measurements ->
automated V_* > 100 km/s cut -> matching (3 controls per counter-rotator,
greedy, no replacement) -> disturbance index -> E = Var(D_CR) - Var(D_CO)
with bootstrap, at k-1, k, k+1 -> registered verdicts.

--dry runs everything but prints only sample counts (no D, no E), for
checking mechanics before the real run.
"""

import argparse
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import manga, manga_measure as mm, models  # noqa: E402
from analysis.g1_manga_select import (MAX_VERR, MIN_HA_SNR, grid_xy,  # noqa: E402
                                      receding_pa)

warnings.filterwarnings("ignore")
SEED = 20261002
N_BOOT = 10_000
V_CUT = 40.0                 # PREREG (g)
PRIMARY = 1.0
ALL_KEYS = ("stellar_vel", "stellar_sigma", "stellar_sigmacorr", "ha_vel",
            "ha_sigma", "ha_instsigma", "ha_flux")


def load_catalogues():
    pa = pd.read_csv(ROOT / "data" / "processed" / "manga_pa.csv")
    pa = pa.drop_duplicates("plateifu", keep="last").set_index("plateifu")
    parent = manga.parent_sample()
    gz = fits.getdata(manga.MANGA / "MaNGA_gz-v2_0_1.fits", 1)
    gz = pd.DataFrame({
        "plateifu": pd.Series(gz["PLATEIFU"].astype(str)).str.strip(),
        "merger": (gz["t08_odd_feature_a24_merger_flag"] == 1)
        | (gz["t08_odd_feature_a21_disturbed_flag"] == 1)}).drop_duplicates(
            "plateifu").set_index("plateifu")
    dl = fits.getdata(manga.MANGA / "manga-morphology-dl-DR17.fits", 1)
    dl = pd.DataFrame({"plateifu": pd.Series(dl["PLATEIFU"].astype(str)).str.strip(),
                       "ttype": dl["T-Type"].astype(float)}).drop_duplicates(
                           "plateifu").set_index("plateifu")
    dap = fits.getdata(manga.MANGA / "dapall-v3_1_1-3.1.0.fits",
                       "HYB10-MILESHC-MASTARSSP")
    dap = pd.DataFrame({"plateifu": pd.Series(dap["PLATEIFU"].astype(str)).str.strip(),
                        "sig1re": dap["STELLAR_SIGMA_1RE"].astype(float)}
                       ).set_index("plateifu")
    d = parent.join(pa, how="inner").join(gz, how="left").join(
        dl, how="left").join(dap, how="left")
    d["merger"] = d["merger"].fillna(False).astype(bool)
    return d


def galaxy_arrays(pifu):
    m = manga.fetch_maps(pifu, keys=ALL_KEYS)
    sv, ss, sc = m["stellar_vel"], m["stellar_sigma"], m["stellar_sigmacorr"]
    hv, hs, hi, hf = m["ha_vel"], m["ha_sigma"], m["ha_instsigma"], m["ha_flux"]
    x, y = grid_xy(sv["value"].shape)

    def err(mp):
        return np.where(mp["ivar"] > 0,
                        1 / np.sqrt(np.clip(mp["ivar"], 1e-30, None)), np.inf)

    es, eg = err(sv), err(hv)
    s_ok = (sv["mask"] == 0) & (sv["ivar"] > 0) & (es <= MAX_VERR)
    g_ok = ((hv["mask"] == 0) & (hv["ivar"] > 0) & (eg <= MAX_VERR)
            & (hf["mask"] == 0)
            & (hf["value"] * np.sqrt(np.clip(hf["ivar"], 0, None)) >= MIN_HA_SNR))
    s_sig2 = ss["value"] ** 2 - sc["value"] ** 2
    g_sig2 = hs["value"] ** 2 - hi["value"] ** 2
    s_sig = np.where((ss["mask"] == 0) & (s_sig2 > 0), np.sqrt(np.abs(s_sig2)), np.nan)
    g_sig = np.where((hs["mask"] == 0) & (g_sig2 > 0), np.sqrt(np.abs(g_sig2)), np.nan)
    # Gas with sigma below instrumental resolution is cold: sigma -> 0.
    g_sig = np.where((hs["mask"] == 0) & (g_sig2 <= 0), 0.0, g_sig)
    nan = np.nan
    maps = dict(vs=np.where(s_ok, sv["value"] - np.nanmedian(sv["value"][s_ok]), nan),
                evs=es, ss=np.where(s_ok, s_sig, nan),
                vg=np.where(g_ok, hv["value"] - np.nanmedian(hv["value"][g_ok]), nan),
                evg=eg, sg=np.where(g_ok, g_sig, nan))
    return maps, x, y, s_ok, g_ok


def measure_one(args):
    pifu, pa_s, pa_g, inc, re_arcsec = args
    try:
        maps, x, y, _, g_ok = galaxy_arrays(pifu)
        rows = mm.measure(maps, x, y, pa_s, pa_g, inc, re_arcsec)
        med_eg = float(np.median(maps["evg"][g_ok])) if g_ok.any() else np.nan
        return pifu, rows, med_eg
    except Exception as e:  # noqa: BLE001
        return pifu, None, repr(e)


def gas_twist(pifu, re_arcsec):
    maps, x, y, _, g_ok = galaxy_arrays(pifu)
    r = np.hypot(x, y)
    out = []
    for sel in (g_ok & (r < re_arcsec), g_ok & (r >= re_arcsec)):
        if sel.sum() < 30:
            return np.nan
        out.append(receding_pa(x[sel], y[sel], maps["vg"][sel],
                               maps["evg"][sel])[0])
    d = abs(out[0] - out[1]) % 360
    return min(d, 360 - d)


def match(cr, pool, cols, n_per=3, seed=SEED):
    allx = pd.concat([cr[cols], pool[cols]])
    mu, sd = allx.mean(), allx.std()
    A = ((cr[cols] - mu) / sd).values
    B = ((pool[cols] - mu) / sd).values
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(cr))
    used = np.zeros(len(pool), bool)
    pairs = {}
    for i in order:
        dist = np.sqrt(((B - A[i]) ** 2).sum(1))
        dist[used] = np.inf
        pick = np.argsort(dist)[:n_per]
        used[pick] = True
        pairs[cr.index[i]] = list(pool.index[pick])
    return pairs


def boot_E(dcr, dco, n=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    E = dcr.var(ddof=1) - dco.var(ddof=1)
    b = np.array([rng.choice(dcr, len(dcr)).var(ddof=1)
                  - rng.choice(dco, len(dco)).var(ddof=1) for _ in range(n)])
    return float(E), float(b.std(ddof=1))


def boot_S(dcr, dco, n=N_BOOT, seed=SEED):
    """S = median(D_CO) - median(D_CR) with bootstrap SD (PREREG g)."""
    rng = np.random.default_rng(seed)
    S = np.median(dco) - np.median(dcr)
    b = np.array([np.median(rng.choice(dco, len(dco)))
                  - np.median(rng.choice(dcr, len(dcr))) for _ in range(n)])
    return float(S), float(b.std(ddof=1))


def s_pred(cr_tab):
    c = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv",
                    index_col=0)
    c = c[c["clean"] & (c["Vflat"] > 0)]
    vals = []
    for _, g in cr_tab.iterrows():
        ref = c[(c["Vflat"] - g["Vc_star"]).abs() < 30]
        if len(ref) == 0:
            continue
        R = g["Re_kpc"]
        vals.append(np.mean(R * models.rt_flow_vorticity(R, ref["omega"].values,
                                                         ref["Rt"].values)))
    return float(np.mean(vals)), len(vals)


def e_pred(cr_tab):
    c = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv",
                    index_col=0)
    c = c[c["clean"] & (c["Vflat"] > 0)]
    vals = []
    for _, g in cr_tab.iterrows():
        ref = c[(c["Vflat"] - g["Vc_star"]).abs() < 30]
        if len(ref) == 0:
            continue
        R = g["Re_kpc"]
        rz = R * models.rt_flow_vorticity(R, ref["omega"].values, ref["Rt"].values)
        vals.append(np.mean(rz ** 2))
    return float(np.mean(vals)), len(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()

    d = load_catalogues()
    ok = (d["e_pa_star"] <= 20) & (d["e_pa_gas"] <= 20) & d["dpa"].notna()
    d = d[ok & ~d["merger"]]
    d["group"] = np.where(d["dpa"] > 150, "CR", np.where(d["dpa"] < 30, "CO", ""))
    d = d[d["group"] != ""].copy()
    d["Re_kpc"] = d["Re_arcsec"] * d["kpc_per_arcsec"]
    log = {"n_CR_candidates": int((d.group == "CR").sum()),
           "n_CO_candidates": int((d.group == "CO").sum())}
    print(f"after PA quality + merger cuts: CR {log['n_CR_candidates']},"
          f" CO {log['n_CO_candidates']}", flush=True)

    jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in d.iterrows()]
    meas, fails = {}, {}
    with ProcessPoolExecutor(args.workers) as ex:
        for pifu, rows, extra in ex.map(measure_one, jobs, chunksize=4):
            if rows is None:
                fails[pifu] = extra
            else:
                meas[pifu] = (rows, extra)
    print(f"measured {len(meas)}, failed {len(fails)}", flush=True)

    rec = []
    for pifu, (rows, med_eg) in meas.items():
        r = {row["r_re"]: row for row in rows}
        base = dict(plateifu=pifu, med_eg=med_eg,
                    gas_resid=np.nanmean([r[c]["gas_resid"] for c in (0.5, 0.75, 1.0, 1.25, 1.5)]))
        for c in (0.5, 1.0, 1.5):
            for k in ("D_lo", "D_mid", "D_hi"):
                base[f"{k}_{c}"] = r[c][k]
        base["Vs_1"] = r[1.0]["Vs"]
        base["Vc_star"] = r[1.0]["Vc_star"]
        rec.append(base)
    t = d.join(pd.DataFrame(rec).set_index("plateifu"), how="inner")
    t = t[np.isfinite(t["D_mid_1.0"]) & (t["Vs_1"] > V_CUT)]
    cols = ["logMstar", "ttype", "log_sig1re", "inc"]
    t["log_sig1re"] = np.log10(t["sig1re"].clip(lower=1))
    t = t.dropna(subset=cols)
    cr, pool = t[t.group == "CR"], t[t.group == "CO"]
    log.update(n_CR_after_cuts=int(len(cr)), n_CO_pool_after_cuts=int(len(pool)))
    print(f"after V_* > {V_CUT:.0f} cut: CR {len(cr)}, CO pool {len(pool)}",
          flush=True)
    if len(cr) == 0 or len(pool) < 3 * len(cr):
        print("insufficient sample; stopping", flush=True)
        return
    pairs = match(cr, pool, cols)
    co = pool.loc[[p for v in pairs.values() for p in v]]
    log["n_CO_matched"] = int(len(co))

    both = pd.concat([cr, co])
    with ProcessPoolExecutor(args.workers) as ex:
        tw = list(ex.map(gas_twist, both.index, both["Re_arcsec"], chunksize=4))
    both["twist"] = tw
    both["resid_ratio"] = both["gas_resid"] / both["med_eg"]
    both["disturb"] = (both["resid_ratio"].rank(pct=True)
                       + both["twist"].rank(pct=True)) / 2
    low = both["disturb"] <= both["disturb"].median()
    print(f"matched: CR {len(cr)}, CO {len(co)};"
          f" low-disturbance half: CR {int(low[cr.index].sum())},"
          f" CO {int(low[co.index].sum())}", flush=True)
    if args.dry:
        print("DRY RUN: mechanics OK; no D or E computed/printed.")
        return

    Ep, n_pred = e_pred(cr)
    Sp, _ = s_pred(cr)
    res = {"log": log, "S_pred": Sp, "E_pred": Ep, "n_pred": n_pred, "radii": {}}
    print(f"\nS_pred = {Sp:.1f} km/s;  E_pred = {Ep:.0f} km^2/s^2 (from {n_pred} CR)")
    for c in (1.0, 0.5, 1.5):
        rr = {}
        for k in ("D_lo", "D_mid", "D_hi"):
            col = f"{k}_{c}"
            a = both.loc[cr.index, col].dropna().values
            b = both.loc[co.index, col].dropna().values
            la = both.loc[cr.index][low[cr.index]][col].dropna().values
            lb = both.loc[co.index][low[co.index]][col].dropna().values
            S, sS = boot_S(a, b)
            Sl, sSl = boot_S(la, lb) if len(la) > 2 and len(lb) > 2 else (np.nan, np.nan)
            E, sE = boot_E(a, b)
            rr[k] = dict(S=S, sS=sS, S_low=Sl, sS_low=sSl, E=E, sE=sE,
                         n_cr=len(a), n_co=len(b), n_cr_low=len(la),
                         med_cr=float(np.median(a)), med_co=float(np.median(b)))
            print(f"  R={c} Re {k:6s}: S = {S:+6.1f} ± {sS:4.1f} km/s"
                  f" (low-dist {Sl:+6.1f} ± {sSl:4.1f}, N_CR={len(la)});"
                  f"  med D_CR {np.median(a):+6.1f}, med D_CO {np.median(b):+6.1f};"
                  f"  E = {E:7.0f} ± {sE:5.0f}")
        res["radii"][str(c)] = rr

    # Registered verdict (PREREG g) at 1 R_e, all three k.
    p = res["radii"]["1.0"]
    if p["D_mid"]["n_cr"] < 40:
        verdict = "UNDERPOWERED (N_CR < 40): inconclusive"
    else:
        fals = all(p[k]["S"] + 2 * p[k]["sS"] < 0.5 * Sp for k in p)
        supp = all(p[k]["S"] > 3 * p[k]["sS"] and abs(p[k]["S"] - Sp) < 2 * p[k]["sS"] for k in p)
        sig_full = all(p[k]["S"] > 3 * p[k]["sS"] for k in p)
        sig_low = all(p[k]["S_low"] > 3 * p[k]["sS_low"] for k in p)
        if fals:
            verdict = "FALSIFIED (stellar-aligned flow at >= 50% of RT-implied amplitude)"
        elif sig_full and not sig_low and p["D_mid"]["n_cr_low"] >= 20:
            verdict = "DISTURBANCE-DRIVEN, not H_vort"
        elif supp:
            verdict = "SUPPORTED"
        else:
            verdict = "INCONCLUSIVE"
    res["verdict"] = verdict
    print(f"\nREGISTERED VERDICT (PREREG g; R = 1 Re, all k): {verdict}")
    both.to_csv(ROOT / "results" / "g1_manga_sample.csv")
    (ROOT / "results" / "g1_manga.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
