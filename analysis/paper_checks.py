#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Numbers quoted in the paper that paper_numbers.py does not produce. POST HOC.

  manga_counts: MaNGA main sample and inclination cut
  c1_counts   : SAMI/CALIFA sample counts (parent split, PA quality, spaxel losses)
  manga       : robust scatter of the slowdowns, counter-rotators without a
                prediction, and the drift term that k multiplies (CR vs controls)
  sami_sigma  : size of the SAMI sample behind the gas-dispersion comparison
  merger      : the zero-twist intercepts with and without the Galaxy Zoo merger
                cut, both computed with the paper's own pipeline (as Table 1)
  seeds       : the intercepts under alternative orders of the greedy matching
  (merger, seeds and pooled quote f and the likelihood ratio against the
   smaller, psi ~ r prediction, P ~ 32 km/s at 1 R_e, i.e. conservatively)
  validation  : synthetic recovery of an injected asymmetry and gas inflow
  pooled      : the zero-twist intercepts with the exploratory (pooled-control)
                definition of the slowdown, for comparison with Table 1
  fs_shear    : exact flowing-space prediction (ODE) against 2 (V_obs - V_bar)
                and against the superseded first-draft formula
  drift_decomposition : medians of the pieces of (sigma_g^2 - sigma_*^2)/(V_g + V_*)
                for counter-rotators and controls ("python paper_checks.py drift"
                recomputes only this block)

The merger check measures the merger-flagged galaxies that the main cache
excludes and stores them in data/processed/manga_meas_mergers.csv.
Output: results/paper_checks.json.
"""
import importlib.util
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.c1_test import MIN_SPAX, PA_OUT, primary  # noqa: E402
from analysis.agm_comparison import model_profiles, predict  # noqa: E402
from analysis.g1_manga_robustness import CACHE, COLS, build  # noqa: E402
from analysis.g1_manga_test import gas_twist, load_catalogues, match, measure_one  # noqa: E402
from src import manga_measure as mm  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MERGER_CACHE = ROOT / "data" / "processed" / "manga_meas_mergers.csv"
PAPER = ROOT / "results" / "paper_manga_cr.csv"
SEEDS = range(1000, 1020)          # alternative greedy-matching orders
N_VALID = 20                       # noise realisations for the synthetic validation


def robust_sd(x):
    x = np.asarray(x, float)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def slowdowns(cr, pool, pairs, radii=(1.0, 1.5), ks=("D_mid",)):
    out = {}
    for c in radii:
        for k in ks:
            col = f"{k}_{c}"
            out[f"s_{k}_{c}"] = pd.Series(
                [np.nanmedian(pool.loc[pairs[n], col]) - cr.loc[n, col]
                 if np.isfinite(pool.loc[pairs[n], col]).any() else np.nan for n in cr.index],
                index=cr.index)
    return out


def intercept(s, twist, pred):
    d = pd.DataFrame({"twist": twist, "s_D_mid_1.5": s, "P": pred})
    d = d[np.isfinite(d["s_D_mid_1.5"]) & np.isfinite(d.twist) & np.isfinite(d.P)]
    return primary(d, float(d.P.mean()), "D_mid")


# ---------------------------------------------------------------- C1 counts
def c1_counts():
    pa = pd.read_csv(PA_OUT)
    good = pa[(pa.e_pa_star <= 20) & (pa.e_pa_gas <= 20) & pa.dpa.notna()]
    sami = pa[pa.survey == "SAMI"]
    return {"parent": int(len(pa)), "parent_by_survey": pa.survey.value_counts().to_dict(),
            "pa_quality": int(len(good)),
            "pa_quality_CR": int((good.dpa > 150).sum()),
            "SAMI_gas_below_min_spaxels": int((sami.n_gas < MIN_SPAX).sum()),
            "SAMI_stars_below_min_spaxels": int((sami.n_star < MIN_SPAX).sum())}


def manga_counts():
    from analysis.g1_manga_select import OUT, parent_sample
    p = parent_sample()
    inc = p[(p["inc"] >= 30) & (p["inc"] <= 80)]
    pa = pd.read_csv(OUT)   # appended on resume, so rows are duplicated; count galaxies
    return {"main_sample": int(len(p)), "inclination_30_80": int(len(inc)),
            "pa_measured": int(pa["plateifu"].drop_duplicates().isin(inc.index).sum()),
            "pa_rows_in_file": int(len(pa))}


# ---------------------------------------------------------------- MaNGA table
def manga(paper):
    out = {}
    for c in ("1.0", "1.5"):
        m = np.isfinite(paper[f"s_D_mid_{c}"]) & np.isfinite(paper.twist)
        nopred = m & ~np.isfinite(paper[f"GM1_{c}"])
        out[c] = {"N_twist": int(m.sum()), "N_fit": int((m & ~nopred).sum()),
                  "robust_sd_s": robust_sd(paper.loc[m, f"s_D_mid_{c}"]),
                  "no_prediction": paper.index[nopred].tolist(),
                  "no_prediction_Vc": paper.Vc_star[nopred].round(1).tolist()}
    # drift term (sigma_g^2 - sigma_*^2)/(V_g + V_*) per unit k = (D_hi - D_lo)/2
    t = build()
    s = t[t.Vs_1 > 40]
    cr, pool = s[s.group == "CR"], s[s.group == "CO"]
    pairs = match(cr, pool, COLS)
    for c in ("1.0", "1.5"):
        x_cr = (cr[f"D_hi_{c}"] - cr[f"D_lo_{c}"]) / 2
        x_co = pd.Series({n: np.nanmedian((pool.loc[pairs[n], f"D_hi_{c}"]
                                           - pool.loc[pairs[n], f"D_lo_{c}"]) / 2) for n in cr.index})
        out[c]["drift_term_median_CR"] = float(np.nanmedian(x_cr))
        out[c]["drift_term_median_controls"] = float(np.nanmedian(x_co))
    return out


def sami_sigma():
    c = pd.read_csv(ROOT / "results" / "c1_sample.csv", index_col=0)
    cr = c[(c.role == "CR") & np.isfinite(c.sg_15)]
    return {"N_CR_with_sigma_1.5": int(len(cr)), "surveys": cr.survey.value_counts().to_dict()}


# ---------------------------------------------------------------- merger cut
def merger_measurements(meas):
    d = load_catalogues()
    d = d[(d.e_pa_star <= 20) & (d.e_pa_gas <= 20) & d.dpa.notna() & d.merger]
    d = d[(d.dpa > 150) | (d.dpa < 30)]
    if MERGER_CACHE.exists():
        extra = pd.read_csv(MERGER_CACHE, index_col=0)
    else:
        jobs = [(p, r.pa_star, r.pa_gas, r.inc, r.Re_arcsec) for p, r in d.iterrows()]
        rec = []
        with ProcessPoolExecutor(8) as ex:
            for pifu, rows, _ in ex.map(measure_one, jobs, chunksize=2):
                if rows is None:
                    continue
                r = {q["r_re"]: q for q in rows}
                o = {"plateifu": pifu, "Vs_1": r[1.0]["Vs"], "Vc_star": r[1.0]["Vc_star"]}
                for c in (1.0, 1.5):
                    for k in ("D_lo", "D_mid", "D_hi"):
                        o[f"{k}_{c}"] = r[c][k]
                rec.append(o)
        extra = pd.DataFrame(rec).set_index("plateifu")
        extra.to_csv(MERGER_CACHE)
    return pd.concat([meas, extra[~extra.index.isin(meas.index)]]), int(len(d))


def manga_table(meas, with_merger_cut):
    """As g1_manga_robustness.build + agm_comparison.manga_table, merger cut optional."""
    d = load_catalogues()
    d = d[(d.e_pa_star <= 20) & (d.e_pa_gas <= 20) & d.dpa.notna()]
    if with_merger_cut:
        d = d[~d.merger]
    d["group"] = np.where(d.dpa > 150, "CR", np.where(d.dpa < 30, "CO", ""))
    d = d[d.group != ""].copy()
    d["Re_kpc"] = d.Re_arcsec * d.kpc_per_arcsec
    d["log_sig1re"] = np.log10(d.sig1re.clip(lower=1))
    t = d.join(meas, how="inner").dropna(subset=COLS)
    t = t[np.isfinite(t["D_mid_1.0"]) & (t.Vs_1 > 40)]
    cr, pool = t[t.group == "CR"].copy(), t[t.group == "CO"]
    for col, ser in slowdowns(cr, pool, match(cr, pool, COLS)).items():
        cr[col] = ser
    with ProcessPoolExecutor(8) as ex:
        cr["twist"] = list(ex.map(gas_twist, cr.index, cr.Re_arcsec))
    return cr, len(pool)


def merger(gm):
    meas = pd.read_csv(CACHE, index_col=0)
    allm, n_flagged = merger_measurements(meas)
    out = {"merger_flagged_CR_or_CO": n_flagged}
    for cut in (True, False):
        cr, npool = manga_table(allm, cut)
        row = {"N_CR": int(len(cr)), "N_CO_pool": int(npool)}
        for c in (1.0, 1.5):
            pred = pd.Series([predict(gm, g.Vc_star, c * g.Re_kpc, "n1") for _, g in cr.iterrows()],
                             index=cr.index)
            row[str(c)] = intercept(cr[f"s_D_mid_{c}"], cr.twist, pred)
        out["with_cut" if cut else "without_cut"] = row
    return out


# ---------------------------------------------------------------- matching seeds
def seeds(paper):
    t = build()
    s = t[t.Vs_1 > 40]
    cr, pool = s[s.group == "CR"], s[s.group == "CO"]
    assert set(cr.index) == set(paper.index), "paper table and build() disagree"
    rows = {"1.0": [], "1.5": []}
    for seed in SEEDS:
        sv = slowdowns(cr, pool, match(cr, pool, COLS, seed=seed))
        for c in ("1.0", "1.5"):
            r = intercept(sv[f"s_D_mid_{c}"], paper.twist, paper[f"GM1_{c}"])
            rows[c].append({k: r[k] for k in ("a", "sd_a", "f_upper95", "log10BF", "N")})
    out = {}
    for c, v in rows.items():
        a = np.array([r["a"] for r in v])
        out[c] = {"seeds": list(SEEDS), "per_seed": v, "a_min": float(a.min()), "a_max": float(a.max()),
                  "a_median": float(np.median(a)),
                  "f_upper95_max": float(max(r["f_upper95"] for r in v)),
                  "log10LR_min": float(min(r["log10BF"] for r in v)),
                  "log10LR_max": float(max(r["log10BF"] for r in v))}
    return out


# ---------------------------------------------------------------- pooled controls
def pooled(paper):
    """Zero-twist intercepts with the exploratory slowdown definition,
    s_i = median(D over ALL matched controls) - D_i, instead of each
    counter-rotator's own three controls (the registered definition used in
    Table 1). Same sample, matching seed, twist and FS prediction."""
    t = build()
    s = t[t.Vs_1 > 40]
    cr, pool = s[s.group == "CR"], s[s.group == "CO"]
    pairs = match(cr, pool, COLS)
    co_all = pool.loc[sorted({p for v in pairs.values() for p in v})]
    out = {"N_CR": int(len(cr)), "N_CO_matched": int(len(co_all))}
    for c in ("1.0", "1.5"):
        sv = co_all[f"D_mid_{c}"].median() - cr[f"D_mid_{c}"]
        r = intercept(sv, paper.twist, paper[f"GM1_{c}"])
        out[c] = {k: r[k] for k in ("a", "sd_a", "f_upper95", "log10BF", "P", "N")}
    return out


# ---------------------------------------------------------------- FS shear correction
def fs_shear(paper):
    """Exact flowing-space (form A) prediction against its approximations.

    From an observed rotation curve V_obs = V_RAR(R) and the baryonic V_bar(R),
    the flow u(R) follows from the prograde root of form A,
        v+ = u - S R/2 + sqrt(S^2 R^2/4 + V_bar^2) = V_obs,   S = u/R - du/dR
    (derivations/01, [A3]). With y = S R/2 this is the ODE
        y = [V_bar^2 - (V_obs - u)^2] / [2 (V_obs - u)],   du/dR = (u - 2 y)/R,
    which is strongly attracting, so the initial condition is immaterial. The
    asymmetry is Delta V = d(R u)/dR = 2 u - 2 y. It is compared with
      (a) 2 (V_obs - V_bar): the first-order-consistent form used in the paper
          (identical to rigid dragging, derivations/01 [A6]);
      (b) d(R u)/dR with u = V_obs - V_bar: the superseded first draft, which
          assumed the additive relation V_obs = V_bar + u that holds only for
          solid-body flow.
    Each is evaluated per counter-rotator with the same SPARC matching as the
    paper's predictions (|V_flat - V_c,*| < 30 km/s, own radius)."""
    from scipy.integrate import solve_ivp
    from src.data import load_mass_models, load_sparc_table, v_bary
    from analysis.s_pred_rar import G_DAGGER
    sp = load_sparc_table()
    sp = sp[(sp.Q < 3) & (sp.Inc >= 30) & (sp.Vflat > 0)]
    prof = {}
    for gid, g in load_mass_models().groupby("ID"):
        if gid not in sp.index:
            continue
        R = g.R.values
        vb = np.clip(v_bary(g.Vgas.values, g.Vdisk.values, g.Vbul.values), 0, None)
        gbar = vb ** 2 / R
        vo = np.sqrt(gbar / (1 - np.exp(-np.sqrt(np.clip(gbar, 1e-6, None) / G_DAGGER))) * R)

        def rhs(r, uu):
            d = max(np.interp(r, R, vo) - uu[0], 1e-6)
            y = (np.interp(r, R, vb) ** 2 - d ** 2) / (2 * d)
            return [(uu[0] - 2 * y) / r]

        sol = solve_ivp(rhs, (R[0], R[-1]), [vo[0] - vb[0]], t_eval=R, method="Radau",
                        rtol=1e-9, atol=1e-9)
        u = sol.y[0]
        d = np.clip(vo - u, 1e-6, None)
        y = (vb ** 2 - d ** 2) / (2 * d)
        vplus = u - y + np.sqrt(y ** 2 + vb ** 2)
        prof[gid] = dict(R=R, exact=2 * u - 2 * y, lin=2 * (vo - vb),
                         old=np.gradient(R * (vo - vb), R), resid=vplus - vo,
                         vflat=sp.loc[gid, "Vflat"])
    out = {"max_abs_root_residual_kms": float(max(np.abs(p["resid"]).max() for p in prof.values()))}
    for c in (1.0, 1.5):
        acc = {k: [] for k in ("exact", "lin", "old")}
        for _, g in paper.iterrows():
            r = c * g.Re_kpc
            sel = [p for p in prof.values()
                   if abs(p["vflat"] - g.Vc_star) < 30 and p["R"].min() <= r <= p["R"].max()]
            if not sel:
                continue
            for k in acc:
                acc[k].append(np.mean([np.interp(r, p["R"], p[k]) for p in sel]))
        acc = {k: np.array(v) for k, v in acc.items()}
        out[str(c)] = {"N": int(len(acc["exact"])),
                       "mean_exact_formA": float(acc["exact"].mean()),
                       "mean_2(Vobs-Vbar)": float(acc["lin"].mean()),
                       "mean_superseded_dRu_dR": float(acc["old"].mean()),
                       "max_abs_exact_minus_2(Vobs-Vbar)_kms": float(np.abs(acc["exact"] - acc["lin"]).max()),
                       "max_abs_exact_minus_superseded_kms": float(np.abs(acc["exact"] - acc["old"]).max())}
    return out


# ---------------------------------------------------------------- drift-term decomposition
def drift_decomposition():
    """Why the shift of D per unit k is more negative in counter-rotators.

    Re-measures the 90 counter-rotators and their 270 matched controls from
    the cached maps and reports group medians at 1.0 and 1.5 R_e of the pieces
    of (sigma_g^2 - sigma_*^2)/(V_g + V_*): the two dispersions, the two
    rotation speeds, the numerator, the denominator and the term itself, plus
    the fitted drift factors. Plain group medians (the paper's -18.5 / -13.8
    km/s are the median over counter-rotators and the median of per-CR control
    medians of (D_hi - D_lo)/2, which include the clipping of k at zero)."""
    t = build()
    s = t[t.Vs_1 > 40]
    cr, pool = s[s.group == "CR"], s[s.group == "CO"]
    pairs = match(cr, pool, COLS)
    co_ids = sorted({p for v in pairs.values() for p in v})
    d = load_catalogues()
    ids = list(cr.index) + co_ids
    jobs = [(p, d.loc[p, "pa_star"], d.loc[p, "pa_gas"], d.loc[p, "inc"], d.loc[p, "Re_arcsec"]) for p in ids]
    rec = {}
    with ProcessPoolExecutor(8) as ex:
        for pifu, rows, _ in ex.map(measure_one, jobs, chunksize=2):
            if rows is not None:
                rec[pifu] = {q["r_re"]: q for q in rows}
    pieces = {
        "sigma_gas": lambda q: q["sg"], "sigma_star": lambda q: q["ss"],
        "V_gas": lambda q: q["Vg"], "V_star": lambda q: q["Vs"],
        "k_gas": lambda q: q["kg"], "k_star": lambda q: q["ks"],
        "numerator_sg2_minus_ss2": lambda q: q["sg"] ** 2 - q["ss"] ** 2,
        "denominator_Vg_plus_Vs": lambda q: q["Vg"] + q["Vs"],
        "term": lambda q: (q["sg"] ** 2 - q["ss"] ** 2) / (q["Vg"] + q["Vs"]),
    }
    out = {"N_CR_measured": int(sum(p in rec for p in cr.index)),
           "N_controls_measured": int(sum(p in rec for p in co_ids))}
    for c in (1.0, 1.5):
        out[str(c)] = {}
        for grp, ids_ in (("CR", list(cr.index)), ("controls", co_ids)):
            vals = {k: [f(rec[p][c]) for p in ids_ if p in rec] for k, f in pieces.items()}
            out[str(c)][grp] = {k: float(np.nanmedian(v)) for k, v in vals.items()}
    return out


# ---------------------------------------------------------------- validation
def validation():
    spec = importlib.util.spec_from_file_location("tmm", ROOT / "tests" / "test_manga_measure.py")
    T = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(T)

    def at_1re(maps, x, y):
        return next(r for r in mm.measure(maps, x, y, T.PA_S, T.PA_G, T.INC, T.RE) if r["r_re"] == 1.0)

    out = {}
    for omega in (3.0, 6.0):
        truth, err = -2 * omega * T.RE, []
        for seed in range(N_VALID):
            T.RNG = np.random.default_rng(seed)
            err.append((at_1re(*T.mock(omega)[:3])["D_mid"] - truth) / abs(truth))
        err = np.array(err)
        out[f"Omega_{omega:g}"] = {"truth": truth, "rel_err_median_abs": float(np.median(np.abs(err))),
                                   "rel_err_max_abs": float(np.abs(err).max())}
    vr, dshift = [], []
    for seed in range(N_VALID):
        T.RNG = np.random.default_rng(seed)
        m0, x, y = T.mock(3.0)
        T.RNG = np.random.default_rng(seed)
        m1, _, _ = T.mock(3.0, vr_gas=-25.0)
        r0, r1 = at_1re(m0, x, y), at_1re(m1, x, y)
        vr.append(r1["VRg"] - 25.0)
        dshift.append(r1["D_mid"] - r0["D_mid"])
    out["inflow_25"] = {"VR_err_max_abs": float(np.abs(vr).max()),
                        "D_shift_max_abs": float(np.abs(dshift).max())}
    return out


def main():
    paper = pd.read_csv(PAPER, index_col=0)
    out_path = ROOT / "results" / "paper_checks.json"
    if sys.argv[1:] == ["drift"]:        # recompute only the drift decomposition
        res = json.loads(out_path.read_text())
        res["drift_decomposition"] = drift_decomposition()
        out_path.write_text(json.dumps(res, indent=2, default=float))
        print_drift(res["drift_decomposition"])
        return
    gm = model_profiles()
    res = {"manga_counts": manga_counts(), "c1_counts": c1_counts(), "manga": manga(paper), "sami_sigma": sami_sigma(),
           "validation": validation(), "seeds": seeds(paper), "merger": merger(gm),
           "pooled": pooled(paper), "fs_shear": fs_shear(paper),
           "drift_decomposition": drift_decomposition()}
    out_path.write_text(json.dumps(res, indent=2, default=float))
    print(f"MaNGA: {res['manga_counts']}")
    c, mg, v, sd, me = res["c1_counts"], res["manga"], res["validation"], res["seeds"], res["merger"]
    print(f"C1: parent {c['parent']} {c['parent_by_survey']}; PA quality {c['pa_quality']}"
          f" ({c['pa_quality_CR']} CR); SAMI below {MIN_SPAX} spaxels: gas"
          f" {c['SAMI_gas_below_min_spaxels']}, stars {c['SAMI_stars_below_min_spaxels']}")
    for k in ("1.0", "1.5"):
        m = mg[k]
        print(f"MaNGA {k} R_e: N(twist) {m['N_twist']}, N(fit) {m['N_fit']}, robust sd(s)"
              f" {m['robust_sd_s']:.1f}; no prediction: V_c = {m['no_prediction_Vc']};"
              f" drift term CR {m['drift_term_median_CR']:+.1f} vs controls"
              f" {m['drift_term_median_controls']:+.1f}")
    print(f"SAMI gas-dispersion sample: {res['sami_sigma']}")
    for k, x in v.items():
        print(f"validation {k}: {x}")
    for k in ("1.0", "1.5"):
        x = sd[k]
        print(f"seeds {k} R_e: a {x['a_min']:+.1f} to {x['a_max']:+.1f} (median {x['a_median']:+.1f});"
              f" max f_up95 {x['f_upper95_max']:.3f}; log10 LR {x['log10LR_min']:+.1f} to {x['log10LR_max']:+.1f}")
    for k in ("with_cut", "without_cut"):
        x = me[k]
        print(f"merger {k}: N_CR {x['N_CR']}, CO pool {x['N_CO_pool']}; "
              + "; ".join(f"{c} R_e a = {x[c]['a']:+.2f} +/- {x[c]['sd_a']:.2f}, f_up95 {x[c]['f_upper95']:.3f},"
                          f" log10 LR {x[c]['log10BF']:+.2f} (N {x[c]['N']})" for c in ("1.0", "1.5")))
    po = res["pooled"]
    print(f"pooled-control slowdown (exploratory definition; N_CR {po['N_CR']}, controls {po['N_CO_matched']}): "
          + "; ".join(f"{c} R_e a = {po[c]['a']:+.2f} +/- {po[c]['sd_a']:.2f}, f_up95 {po[c]['f_upper95']:.3f},"
                      f" log10 LR {po[c]['log10BF']:+.2f} (N {po[c]['N']})" for c in ("1.0", "1.5")))
    fsh = res["fs_shear"]
    for c in ("1.0", "1.5"):
        x = fsh[c]
        print(f"flowing space {c} R_e: exact {x['mean_exact_formA']:.2f}; 2(Vobs-Vbar) {x['mean_2(Vobs-Vbar)']:.2f}"
              f" (max |diff| per CR {x['max_abs_exact_minus_2(Vobs-Vbar)_kms']:.2f} km/s);"
              f" superseded d(Ru)/dR {x['mean_superseded_dRu_dR']:.2f}"
              f" (max |diff| {x['max_abs_exact_minus_superseded_kms']:.2f} km/s); N {x['N']};"
              f" root residual {fsh['max_abs_root_residual_kms']:.1e}")
    print_drift(res["drift_decomposition"])

def print_drift(dd):
    for c in ("1.0", "1.5"):
        for grp in ("CR", "controls"):
            x = dd[c][grp]
            print(f"drift pieces {c} R_e {grp:8s}: sigma_g {x['sigma_gas']:.1f}, sigma_* {x['sigma_star']:.1f},"
                  f" V_g {x['V_gas']:.1f}, V_* {x['V_star']:.1f}, k_g {x['k_gas']:.2f}, k_* {x['k_star']:.2f};"
                  f" numerator {x['numerator_sg2_minus_ss2']:+.0f}, denominator {x['denominator_Vg_plus_Vs']:.0f},"
                  f" term {x['term']:+.1f}")


if __name__ == "__main__":
    main()
