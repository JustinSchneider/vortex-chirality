#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Power of the CONFIRMATORY design C1 (SAMI DR3 + CALIFA DR3), calibrated on
the MaNGA exploratory counter-rotators (results/g1_manga_discriminant_cr.csv).

Per counter-rotator i at 1.5 R_e: slowdown s_i = median(D_CO, matched) - D_CR,i,
gas kinematic twist t_i (deg). Model s_i = a + b t_i + e_i, fitted by
Theil-Sen; a is the slowdown of a counter-rotator with settled gas.
  H_C (conventional): a = 0
  H_V (vortex):       a = P = mean_i R*zeta_i (RAR normalisation, 1.5 R_e)
Evidence: BF_VC = exp[(a^2 - (a - P)^2) / (2 sigma_a^2)], sigma_a the
bootstrap SD of the intercept. Simulations resample the MaNGA twist
distribution and Theil-Sen residuals; the slope b is held at the MaNGA value
under both hypotheses.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEED = 20261002
N_BOOT = 400
N_SIM = 1000


def theil_sen(x, y):
    """Vectorised Theil-Sen over the leading axis: x, y shape (B, N)."""
    x, y = np.atleast_2d(x), np.atleast_2d(y)
    i, j = np.triu_indices(x.shape[1], 1)
    dx, dy = x[:, j] - x[:, i], y[:, j] - y[:, i]
    with np.errstate(invalid="ignore", divide="ignore"):
        sl = np.where(dx != 0, dy / dx, np.nan)
    b = np.nanmedian(sl, axis=1)
    a = np.median(y - b[:, None] * x, axis=1)
    return a, b


def intercept_and_sd(x, y, rng, n_boot=N_BOOT):
    a, b = theil_sen(x, y)
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    ab, _ = theil_sen(x[idx], y[idx])
    return float(a[0]), float(b[0]), float(np.nanstd(ab))


def log10_bf(a, sd, P):
    return (a ** 2 - (a - P) ** 2) / (2 * sd ** 2) / np.log(10)


def main():
    c = pd.read_csv(ROOT / "results" / "g1_manga_discriminant_cr.csv", index_col=0)
    m = c[["s_15", "twist", "rz_15"]].dropna()
    x, y = m.twist.values.astype(float), m.s_15.values
    rng = np.random.default_rng(SEED)
    a0, b0, sd0 = intercept_and_sd(x, y, rng)
    P = float(m.rz_15.mean())
    res = y - (a0 + b0 * x)
    res = res - np.median(res)
    print(f"MaNGA calibration (N = {len(m)}): intercept {a0:+.1f} +/- {sd0:.1f} km/s,"
          f" slope {b0:.2f} km/s/deg, robust residual SD {1.4826 * np.median(np.abs(res)):.0f},"
          f" P(vortex) = {P:.1f}; log10 BF_VC = {log10_bf(a0, sd0, P):+.2f}")
    out = {"calibration": dict(N=len(m), a=a0, sd_a=sd0, b=b0, P=P,
                               log10BF=float(log10_bf(a0, sd0, P)))}
    print(f"{'N':>4} {'truth':>5}  P(BF>10 right) P(3<BF<10 right) P(BF>10 wrong) P(BF>3 wrong)")
    for N in (10, 15, 20, 25, 30):
        for truth, A in (("C", 0.0), ("V", P)):
            lb = []
            for _ in range(N_SIM):
                xs = rng.choice(x, N)
                ys = A + b0 * xs + rng.choice(res, N)
                a, _, sd = intercept_and_sd(xs, ys, rng)
                lb.append(log10_bf(a, sd, P))
            lb = np.array(lb) * (1 if truth == "V" else -1)  # >0 favours truth
            r = dict(dec=float(np.mean(lb > 1)), mod=float(np.mean((lb > 0.5) & (lb <= 1))),
                     wrong_dec=float(np.mean(lb < -1)), wrong_mod=float(np.mean(lb < -0.5)))
            out[f"{N}_{truth}"] = r
            print(f"{N:>4} {truth:>5}  {r['dec']:14.2f} {r['mod']:16.2f} {r['wrong_dec']:14.3f}"
                  f" {r['wrong_mod']:13.3f}", flush=True)
    (ROOT / "results" / "c1_power.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
