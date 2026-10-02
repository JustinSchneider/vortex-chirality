#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Power check for the MaNGA Gate G1 (PREREG amendment c), on synthetic data.

The per-galaxy statistic D = gas minus drift-corrected stellar circular
speed at 1 R_e. The simulation models:
  - signal: counter-rotators (CR) get D = ±R*zeta (random sign), with R*zeta
    drawn from clean SPARC galaxies (V_flat > 100) evaluated at R = R_eff;
  - nuisance shared by both groups: gas disequilibrium/non-circular scatter
    (ATLAS3D CO vs JAM: 15% RMS), per-galaxy k error (sigma_k = 0.5), and
    measurement error;
  - extra disturbance in CR gas (factor f_dist on the nuisance), the main
    false-positive risk for a variance test.
Primary statistic E = Var(D_CR) - Var(D_CO), bootstrap uncertainty.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src import models  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RNG = np.random.default_rng(20261002)
N_SIM = 400
N_BOOT = 400
V_TYP, SIG_TYP = 180.0, 80.0          # typical MaNGA fast rotator at 1 R_e
GAS_SCATTER_FRAC = 0.15               # Davis+2013 CO vs JAM RMS
SIG_K = 0.5
MEAS = 7.0                            # km/s, ring-averaged velocity error


def expected_rzeta():
    c = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv",
                    index_col=0)
    c = c[c["clean"] & (c["Vflat"] > 100)]
    return R_zeta_at(c["Reff"].values, c["omega"].values, c["Rt"].values)


def R_zeta_at(R, om, Rt):
    return R * models.rt_flow_vorticity(R, om, Rt)


def nuisance(n, f=1.0):
    gas = RNG.normal(0, GAS_SCATTER_FRAC * V_TYP, n)
    # k error moves the stellar drift term k sigma^2 / (2V)
    kerr = RNG.normal(0, SIG_K, n) * SIG_TYP**2 / (2 * V_TYP)
    meas = RNG.normal(0, MEAS * np.sqrt(2), n)
    return f * gas + kerr + meas


def one_trial(n_cr, rz_pool, signal=True, f_dist=1.0):
    n_co = 3 * n_cr
    d_co = nuisance(n_co)
    rz = RNG.choice(rz_pool, n_cr) if signal else np.zeros(n_cr)
    d_cr = RNG.choice([-1, 1], n_cr) * rz + nuisance(n_cr, f_dist)
    E = d_cr.var(ddof=1) - d_co.var(ddof=1)
    boots = np.array([
        RNG.choice(d_cr, n_cr).var(ddof=1) - RNG.choice(d_co, n_co).var(ddof=1)
        for _ in range(N_BOOT)])
    return E, boots.std(ddof=1)


def main():
    pool = expected_rzeta()
    E_pred = float(np.mean(pool**2))
    print("=" * 72)
    print("MaNGA G1 power check (synthetic data only)")
    print("=" * 72)
    print(f"Expected R*zeta at R_eff (SPARC, V_flat > 100, N = {len(pool)}):"
          f" median {np.median(pool):.0f}, RMS {np.sqrt(E_pred):.0f} km/s")
    print(f"Predicted E = <(R zeta)^2> = {E_pred:.0f} km^2/s^2")
    print(f"Shared nuisance RMS per galaxy ≈"
          f" {nuisance(200000).std():.0f} km/s\n")

    res = {"E_pred": E_pred, "rows": []}
    print(f"{'N_CR':>5} {'f_dist':>6} | {'detect 3σ (signal)':>19}"
          f" | {'falsify (no signal)':>20} | {'false +3σ (no signal)':>22}")
    for n_cr in (20, 40, 60, 100, 150):
        for f in (1.0, 1.3):
            det = fal = fp = 0
            for _ in range(N_SIM):
                E, s = one_trial(n_cr, pool, True, f)
                det += E > 3 * s
                E0, s0 = one_trial(n_cr, pool, False, f)
                fal += (E0 + 2 * s0) < 0.25 * E_pred
                fp += E0 > 3 * s0
            row = dict(n_cr=n_cr, f_dist=f, p_detect=det / N_SIM,
                       p_falsify_null=fal / N_SIM, p_false_pos=fp / N_SIM)
            res["rows"].append(row)
            print(f"{n_cr:5d} {f:6.1f} | {det / N_SIM:19.2f} |"
                  f" {fal / N_SIM:20.2f} | {fp / N_SIM:22.2f}")
    (ROOT / "results" / "g1_manga_power.json").write_text(
        json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
