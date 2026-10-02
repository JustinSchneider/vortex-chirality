#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""MaNGA G1, stage 1: kinematic position angles and gas-star misalignment.

Implements PREREG amendments (c) and (e). Fetches stellar velocity, Halpha
velocity and Halpha flux maps for the parent sample, fits receding-side
kinematic PAs, and writes data/processed/manga_pa.csv with PAs, errors,
spaxel counts and Delta PA only (no rotation amplitudes).

Resumable: maps are cached, and galaxies already in the output are skipped.
Usage: python analysis/g1_manga_select.py [--limit N] [--workers W]
"""

import argparse
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.manga import SPAXEL_ARCSEC, fetch_maps, parent_sample  # noqa: E402

warnings.filterwarnings("ignore")
OUT = ROOT / "data" / "processed" / "manga_pa.csv"
MIN_SPAX = 50          # counted after decimation
MAX_VERR = 30.0
MIN_HA_SNR = 5.0
# PA fits use every 2nd spaxel in x and y: adjacent spaxels are correlated
# (PSF FWHM ~2.5" = 5 spaxels), and this cuts fit time ~4x (PREREG e note).
DECIMATE = 2


def receding_pa(x, y, v, verr):
    """Receding-side kinematic PA in [0, 360) and its error (degrees)."""
    from pafit.fit_kinematic_pa import fit_kinematic_pa
    v = v - np.median(v)
    ang, err, _ = fit_kinematic_pa(x, y, v, quiet=True, plot=False,
                                   dvel=np.clip(verr, 5, None), nsteps=181)
    # pafit measures PA counter-clockwise from +y (towards -x); the pipeline
    # (src/manga_measure.disk_coords) measures from +y towards +x.
    ang = (-ang) % 180.0
    th = np.radians(ang)
    proj = x * np.sin(th) + y * np.cos(th)   # distance along the PA axis
    if np.sum(v * proj) < 0:
        ang = (ang + 180.0) % 360.0
    return float(ang), float(err)


def grid_xy(shape):
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    return ((xx - (nx - 1) / 2) * SPAXEL_ARCSEC,
            (yy - (ny - 1) / 2) * SPAXEL_ARCSEC)


_SESSION = None


def process(pifu):
    global _SESSION
    if _SESSION is None:
        _SESSION = requests.Session()
    m = fetch_maps(pifu, session=_SESSION)
    sv, hv, hf = m["stellar_vel"], m["ha_vel"], m["ha_flux"]
    x, y = grid_xy(sv["value"].shape)
    s_err = np.where(sv["ivar"] > 0, 1 / np.sqrt(np.clip(sv["ivar"], 1e-30, None)), np.inf)
    g_err = np.where(hv["ivar"] > 0, 1 / np.sqrt(np.clip(hv["ivar"], 1e-30, None)), np.inf)
    g_snr = hf["value"] * np.sqrt(np.clip(hf["ivar"], 0, None))
    s_ok = (sv["mask"] == 0) & (sv["ivar"] > 0) & (s_err <= MAX_VERR)
    g_ok = ((hv["mask"] == 0) & (hv["ivar"] > 0) & (g_err <= MAX_VERR)
            & (hf["mask"] == 0) & (g_snr >= MIN_HA_SNR))
    dec = np.zeros_like(s_ok)
    dec[::DECIMATE, ::DECIMATE] = True
    s_ok &= dec
    g_ok &= dec
    row = {"plateifu": pifu, "n_star": int(s_ok.sum()), "n_gas": int(g_ok.sum())}
    if row["n_star"] >= MIN_SPAX:
        row["pa_star"], row["e_pa_star"] = receding_pa(
            x[s_ok], y[s_ok], sv["value"][s_ok], s_err[s_ok])
    if row["n_gas"] >= MIN_SPAX:
        row["pa_gas"], row["e_pa_gas"] = receding_pa(
            x[g_ok], y[g_ok], hv["value"][g_ok], g_err[g_ok])
    if "pa_star" in row and "pa_gas" in row:
        d = abs(row["pa_gas"] - row["pa_star"]) % 360.0
        row["dpa"] = min(d, 360.0 - d)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    parent = parent_sample()
    parent = parent[(parent["inc"] >= 30) & (parent["inc"] <= 80)]
    done = set(pd.read_csv(OUT)["plateifu"]) if OUT.exists() else set()
    todo = [p for p in parent.index if p not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"parent {len(parent)}, done {len(done)}, to process {len(todo)}",
          flush=True)

    rows, n_fail = [], 0
    with ProcessPoolExecutor(args.workers) as ex:
        futs = {ex.submit(process, p): p for p in todo}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                rows.append(f.result())
            except Exception as e:  # network or map error: retried next run
                n_fail += 1
                print(f"  fail {futs[f]}: {type(e).__name__}", flush=True)
            if i % 100 == 0 or i == len(todo):
                pd.DataFrame(rows).to_csv(OUT, mode="a", index=False,
                                          header=not OUT.exists())
                rows = []
                print(f"  {i}/{len(todo)} processed ({n_fail} failed)",
                      flush=True)


if __name__ == "__main__":
    main()
