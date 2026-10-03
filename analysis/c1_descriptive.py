#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""C1 secondary quantities, DESCRIPTIVE ONLY (N < 10; see PREREG.md, 2026-10-02).
Uses the frozen functions of analysis/c1_test.py on results/c1_sample.csv."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.c1_test import balance, boot_S, primary  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

t = pd.read_csv(ROOT / "results" / "c1_sample.csv", index_col=0)
cr, co = t[t.role == "CR"], t[t.role == "CO"]
prim = cr[np.isfinite(cr["s_D_mid_1.5"]) & np.isfinite(cr.twist)]
P = float(np.nanmean(prim.rz_15))
out = {"N_primary": int(len(prim)), "P": P,
       "per_CR": cr[["survey", "twist", "s_D_lo_1.0", "s_D_mid_1.0", "s_D_hi_1.0",
                     "s_D_lo_1.5", "s_D_mid_1.5", "s_D_hi_1.5", "rz_1", "rz_15", "Vc_star"]]
       .round(1).reset_index().to_dict("records")}
out["primary_descriptive"] = {k: primary(prim, P, k) for k in ("D_mid", "D_lo", "D_hi")}
for c in (1.0, 1.5):
    for k in ("D_lo", "D_mid", "D_hi"):
        a, b = cr[f"{k}_{c}"].dropna().values, co[f"{k}_{c}"].dropna().values
        if len(a) > 2 and len(b) > 2:
            S, sS = boot_S(a, b)
            out[f"S_{k}_{c}"] = dict(S=S, sS=sS, n_cr=len(a), n_co=len(b))
out["balance"] = balance(cr, co)
print(json.dumps(out, indent=1, default=float))
(ROOT / "results" / "c1_descriptive.json").write_text(json.dumps(out, indent=2, default=float))
