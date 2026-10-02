#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""S2c (PREREG h): does R_t/R_d depend on M_BH/M_bar at fixed Pi2?
Dynamical M_BH from van den Bosch 2016 (Davis+2017 adds no new SPARC match).
Upper limits (lower uncertainty >= 1 dex) excluded. N < 15 -> EXPLORATORY."""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261002)

def norm(n):
    n = n.upper().replace(" ", "")
    m = re.match(r"^(NGC|UGC|IC|ESO|DDO|F|PGC)0*(\d.*)$", n)
    return m.group(1) + m.group(2) if m else n

rows = []
for f in ("table2.dat", "table3.dat"):
    for l in open(ROOT / "data" / "bh" / f"vdb16_{f}"):
        rows.append(dict(Name=l[0:14].strip(), Dist=float(l[17:25]), logM=float(l[31:36]),
                         e=float(l[42:46]), Meth=l[104:110].strip()))
bh = pd.DataFrame(rows)
bh["key"] = bh.Name.map(norm)
bh = bh[bh.e < 1.0]                                   # exclude upper limits
c = pd.read_csv(ROOT / "data" / "processed" / "sparc_canonical.csv", index_col=0)
c["key"] = [norm(n) for n in c.index]
m = c.reset_index().merge(bh, on="key").set_index("Galaxy")
m["logMBH"] = m.logM + np.log10(m.D / m.Dist)          # rescale to SPARC distance
m["Pi3"] = 10 ** m.logMBH / m.M_bar

def partial(df, y):
    r = lambda s: stats.rankdata(s)
    zx, zy, zz = r(np.log10(df.Pi3)), r(df[y]), r(df.Pi2)
    res = lambda a: a - np.polyval(np.polyfit(zz, a, 1), zz)
    ex, ey = res(zx), res(zy)
    rho = stats.pearsonr(ex, ey).statistic
    null = [stats.pearsonr(ex, RNG.permutation(ey)).statistic for _ in range(10_000)]
    return float(rho), float((np.sum(np.abs(null) >= abs(rho)) + 1) / 10_001)

out = {}
for name, s in (("clean", m[m.clean]), ("taper_pref", m[m.taper_pref]), ("all_ok", m[m.ok_fit])):
    print(f"\n[{name}] N = {len(s)}: {', '.join(s.index)}")
    out[name] = {"N": len(s), "galaxies": list(s.index)}
    if len(s) >= 4:
        for y in ("Rt_over_Rd", "Rt", "omega"):
            rho, p = partial(s, y)
            out[name][y] = dict(rho=rho, p=p)
            print(f"   partial rho({y}, M_BH/M_bar | Pi2) = {rho:+.2f}, p = {p:.2f}")
        raw = stats.spearmanr(s.logMBH, s.Rt_over_Rd)
        print(f"   raw Spearman(Rt/Rd, log M_BH) = {raw.statistic:+.2f}, p = {raw.pvalue:.2f}")
verdict = "EXPLORATORY (N < 15): no test of H_seed possible" if out["clean"]["N"] < 15 else "see stats"
print(f"\nS2c status: {verdict}")
out["verdict"] = verdict
(ROOT / "results" / "s2c_seed.json").write_text(json.dumps(out, indent=2))
