#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Numbers and figures for the paper, from one consistent MaNGA table
(results/agm_comparison_manga_cr.csv: the 90 registered counter-rotators,
per-CR matched controls) plus the C1 registered sample. POST HOC except where
it restates registered results.

Models (all normalised to the RAR on matched SPARC galaxies, at each
counter-rotator's own radius):
  GM1: strong-gravitomagnetic GR, psi = C0 r:  Delta V = (V^2 - V_bar^2)/V
  GM2: flowing space (any flow profile, to first order in the shear), the
       vector-potential form, and rigid dragging psi = K r^2:
       Delta V = 2 (V - V_bar)                    (derivations/01 [A6], 04)
  FS_superseded: the first draft's flowing-space formula, Delta V = d(R u)/dR
       with u = V_RAR - V_bar. It assumed V_obs = V_bar + u, which holds only
       for solid-body flow (wrong at first order in the shear). Kept in the
       JSON and the CSV for the record; not used in the paper. The C1
       registration's P was computed with it (analysis/c1_test.py, rz_15).
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.agm_comparison import model_profiles, predict  # noqa: E402
from analysis.c1_power import theil_sen  # noqa: E402
from analysis.c1_test import primary  # noqa: E402
from analysis.g1_manga_discriminant_exploratory import pred_rz  # noqa: E402
from analysis.s_pred_rar import rar_profiles  # noqa: E402

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
FIG = ROOT / "manuscript" / "figures"
RNG = np.random.default_rng(20261003)
MODELS = {"GM1": r"GR dragging, $\psi\propto r$",
          "GM2": r"flowing space / rigid dragging, $\psi\propto r^2$"}
SUPERSEDED = "FS"      # columns FS_1.0 / FS_1.5 are kept for the record only


def add_predictions(t, gm, fs):
    for c, tag in ((1.0, "1.0"), (1.5, "1.5")):
        t[f"FS_{tag}"] = [pred_rz(fs, g.Vc_star, c * g.Re_kpc) for _, g in t.iterrows()]
        t[f"GM1_{tag}"] = [predict(gm, g.Vc_star, c * g.Re_kpc, "n1") for _, g in t.iterrows()]
        t[f"GM2_{tag}"] = [predict(gm, g.Vc_star, c * g.Re_kpc, "n2") for _, g in t.iterrows()]
    return t


def boot_median(x, n=5000):
    x = np.asarray(x)
    b = np.median(RNG.choice(x, (n, len(x))), axis=1)
    return float(np.median(x)), float(b.std())


def main():
    gm, fs = model_profiles(), rar_profiles()
    t = pd.read_csv(ROOT / "results" / "agm_comparison_manga_cr.csv", index_col=0)
    t = add_predictions(t, gm, fs)
    out = {"N_CR": int(len(t))}
    # correlations and twist split
    for tag in ("1.0", "1.5"):
        s = t[f"s_D_mid_{tag}"]
        m = np.isfinite(s) & np.isfinite(t.twist)
        r = stats.spearmanr(t.twist[m], s[m])
        mp = m & np.isfinite(t[f"GM1_{tag}"])
        rp = stats.spearmanr(t[f"GM1_{tag}"][mp], s[mp])
        med = t.twist[m].median()
        lo, hi = s[m & (t.twist <= med)], s[m & (t.twist > med)]
        out[f"R{tag}"] = dict(
            N=int(m.sum()), rho_twist=float(r.statistic), p_twist=float(r.pvalue),
            rho_pred=float(rp.statistic), p_pred=float(rp.pvalue),
            twist_median=float(med), low_twist=boot_median(lo), high_twist=boot_median(hi),
            n_low=int(len(lo)), n_high=int(len(hi)),
            pred={k: float(t[f"{k}_{tag}"][m].mean()) for k in MODELS},
            pred_FS_superseded=float(t[f"{SUPERSEDED}_{tag}"][m].mean()))
        out[f"R{tag}"]["intercepts"] = {}
        for key in MODELS:
            for k in ("D_lo", "D_mid", "D_hi"):
                d = pd.DataFrame({"twist": t.twist, f"s_{k}_1.5": t[f"s_{k}_{tag}"],
                                  "P": t[f"{key}_{tag}"]})
                d = d[np.isfinite(d[f"s_{k}_1.5"]) & np.isfinite(d.twist) & np.isfinite(d.P)]
                out[f"R{tag}"]["intercepts"][f"{key}_{k}"] = primary(d, float(d.P.mean()), k)
        print(f"R = {tag} R_e: N = {m.sum()}, rho(s, twist) = {r.statistic:+.2f} (p = {r.pvalue:.3f}),"
              f" rho(s, GM1 pred) = {rp.statistic:+.2f} (p = {rp.pvalue:.2f});"
              f" low-twist median s = {out[f'R{tag}']['low_twist'][0]:+.1f} +/- {out[f'R{tag}']['low_twist'][1]:.1f},"
              f" high-twist {out[f'R{tag}']['high_twist'][0]:+.1f} +/- {out[f'R{tag}']['high_twist'][1]:.1f}")
        print(f"   (superseded first-draft FS formula would give P = {out[f'R{tag}']['pred_FS_superseded']:.1f})")
        for key in MODELS:
            row = [out[f"R{tag}"]["intercepts"][f"{key}_{k}"] for k in ("D_lo", "D_mid", "D_hi")]
            print(f"   {key}: P = {row[1]['P']:.1f}; a = {row[1]['a']:+.1f} +/- {row[1]['sd_a']:.1f};"
                  f" f_up95 = {row[1]['f_upper95']:.2f}; log10 LR (k-1,k,k+1) = "
                  + ", ".join(f"{x['log10BF']:+.1f}" for x in row))
    # SAMI registered sample (descriptive)
    c1 = pd.read_csv(ROOT / "results" / "c1_sample.csv", index_col=0)
    c1 = c1[(c1.role == "CR") & np.isfinite(c1["s_D_mid_1.5"]) & np.isfinite(c1.twist)]
    c1["GM1_1.5"] = [predict(gm, g.Vc_star, 1.5 * g.Re_kpc, "n1") for _, g in c1.iterrows()]
    c1["GM2_1.5"] = [predict(gm, g.Vc_star, 1.5 * g.Re_kpc, "n2") for _, g in c1.iterrows()]
    out["SAMI"] = dict(twist=c1.twist.tolist(), s=c1["s_D_mid_1.5"].tolist(),
                       GM1=float(c1["GM1_1.5"].mean()), GM2=float(c1["GM2_1.5"].mean()),
                       P_registered_FS_superseded=float(c1.rz_15.mean()))
    # radial profile of predictions for the median counter-rotator
    vc, re = float(t.Vc_star.median()), float(t.Re_kpc.median())
    radii = np.linspace(0.5, 4.0, 15)
    prof = {k: [] for k in MODELS}
    for m in radii:
        prof["GM1"].append(predict(gm, vc, m * re, "n1"))
        prof["GM2"].append(predict(gm, vc, m * re, "n2"))
    out["profile"] = dict(Vc=vc, Re_kpc=re, radii=radii.tolist(), **{k: list(map(float, v)) for k, v in prof.items()})
    (ROOT / "results" / "paper_numbers.json").write_text(json.dumps(out, indent=2, default=float))
    t.to_csv(ROOT / "results" / "paper_manga_cr.csv")
    figures(t, out, c1)


def figures(t, out, c1):
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})
    col = {"GM1": "#d95f02", "GM2": "#7570b3"}
    # Fig. 1: predicted asymmetry against radius
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    p = out["profile"]
    for k, lab in MODELS.items():
        ax.plot(p["radii"], p[k], color=col[k], lw=1.4, label=lab)
    for n, (tag, x) in enumerate((("1.0", 1.0), ("1.5", 1.5))):
        ic = out[f"R{tag}"]["intercepts"]
        a = ic["GM1_D_mid"]
        # drift-factor systematic: range of the intercept over k-1, k, k+1
        lo = min(ic[f"GM1_{k}"]["a"] for k in ("D_lo", "D_mid", "D_hi"))
        hi = max(ic[f"GM1_{k}"]["a"] for k in ("D_lo", "D_mid", "D_hi"))
        ax.fill_between([x - 0.09, x + 0.09], lo, hi, color="0.85", lw=0,
                        label=r"drift-factor range ($k\pm1$)" if n == 0 else None)
        ax.errorbar(x, a["a"], yerr=1.96 * a["sd_a"], fmt="none", ecolor="k", lw=0.7, capsize=2,
                    label="95% interval" if n == 0 else None)
        ax.errorbar(x, a["a"], yerr=a["sd_a"], fmt="o", color="k", ms=4, lw=1.8, capsize=0,
                    label=r"MaNGA zero-twist slowdown ($\pm1\sigma$)" if n == 0 else None)
    ax.axhline(0, color="0.5", lw=0.6, ls=":")
    ax.text(3.95, 2, "no velocity-dependent force", ha="right", va="bottom", fontsize=6, color="0.4")
    ax.set_xlabel(r"$R/R_{\rm e}$")
    ax.set_ylabel(r"$\Delta V = v_{\rm pro}-|v_{\rm retro}|$ (km s$^{-1}$)")
    ax.set_xlim(0.4, 4.1)
    ax.set_ylim(-45, 100)
    h, lab = ax.get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, fontsize=6, loc="lower center", ncol=2,
               columnspacing=1.0, handlelength=1.6, handletextpad=0.5)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    fig.savefig(FIG / "fig_predictions.pdf")
    plt.close(fig)
    # Fig. 2: slowdown against twist, with binned medians; sqrt twist axis
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.3), sharey=True)
    for ax, tag in zip(axes, ("1.0", "1.5")):
        m = np.isfinite(t[f"s_D_mid_{tag}"]) & np.isfinite(t.twist)
        x, y = t.twist[m].values.astype(float), t[f"s_D_mid_{tag}"][m].values
        ax.set_xscale("function", functions=(np.sqrt, np.square))
        ax.scatter(x, y, s=8, c="0.65", lw=0, label="individual MaNGA galaxies")
        edges = np.quantile(x, [0, 0.25, 0.5, 0.75, 1.0])
        bx, by, be = [], [], []
        for q in range(4):
            sel = (x >= edges[q]) & ((x < edges[q + 1]) if q < 3 else (x <= edges[q + 1]))
            if sel.sum() >= 5:
                med, err = boot_median(y[sel])
                bx.append(float(np.median(x[sel])))
                by.append(med)
                be.append(err)
        ax.errorbar(bx, by, yerr=be, fmt="s", color="k", ms=4.5, capsize=2, lw=1,
                    label="median of each quarter of the sample")
        a, b = theil_sen(x, y)
        xx = np.linspace(0, 180, 200)
        ax.plot(xx, a[0] + b[0] * xx, color="k", lw=0.6, ls="-", alpha=0.6, label="Theil-Sen fit")
        for k in ("GM1", "GM2"):
            ax.axhline(out[f"R{tag}"]["pred"][k], color=col[k], lw=1.1, ls="--",
                       label=MODELS[k] + " (prediction)")
        if tag == "1.5":
            ax.scatter(out["SAMI"]["twist"], out["SAMI"]["s"], s=24, marker="D",
                       facecolor="none", edgecolor="#e7298a", lw=0.9, label="SAMI (registered, N=5)")
        ax.axhline(0, color="0.4", lw=0.7, ls=":")
        ax.set_xticks([0, 5, 10, 25, 50, 100, 180])
        ax.set_xticklabels(["0", "5", "10", "25", "50", "100", "180"])
        ax.set_xlim(0, 185)
        ax.set_xlabel("gas kinematic twist (deg; square-root scale)")
        ax.set_title(rf"$R = {tag}\,R_{{\rm e}}$", fontsize=8)
        ax.set_ylim(-230, 270)
    axes[0].set_ylabel(r"slowdown $s$ (km s$^{-1}$)")
    h, lab = axes[1].get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, fontsize=6.5, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(FIG / "fig_twist.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
