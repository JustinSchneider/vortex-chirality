#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Numbers and figures for the paper, from one consistent MaNGA table
(results/agm_comparison_manga_cr.csv: the 90 registered counter-rotators,
per-CR matched controls) plus the C1 registered sample. POST HOC except where
it restates registered results.

Models (all normalised to the RAR on matched SPARC galaxies, at each
counter-rotator's own radius):
  FS : flowing space / vector potential, Delta V = R zeta = d(R u)/dR, u = V_RAR - V_bar
  GM1: strong-gravitomagnetic GR, psi = C0 r:  Delta V = (V^2 - V_bar^2)/V
  GM2: strong-gravitomagnetic GR, psi = K r^2: Delta V = 2 (V - V_bar)
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
MODELS = {"FS": "flowing space", "GM1": r"GR dragging, $\psi\propto r$",
          "GM2": r"GR dragging, $\psi\propto r^2$"}


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
            pred={k: float(t[f"{k}_{tag}"][m].mean()) for k in MODELS})
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
        for key in MODELS:
            row = [out[f"R{tag}"]["intercepts"][f"{key}_{k}"] for k in ("D_lo", "D_mid", "D_hi")]
            print(f"   {key}: P = {row[1]['P']:.1f}; a = {row[1]['a']:+.1f} +/- {row[1]['sd_a']:.1f};"
                  f" f_up95 = {row[1]['f_upper95']:.2f}; log10 LR (k-1,k,k+1) = "
                  + ", ".join(f"{x['log10BF']:+.1f}" for x in row))
    # SAMI registered sample (descriptive)
    c1 = pd.read_csv(ROOT / "results" / "c1_sample.csv", index_col=0)
    c1 = c1[(c1.role == "CR") & np.isfinite(c1["s_D_mid_1.5"]) & np.isfinite(c1.twist)]
    c1["GM1_1.5"] = [predict(gm, g.Vc_star, 1.5 * g.Re_kpc, "n1") for _, g in c1.iterrows()]
    out["SAMI"] = dict(twist=c1.twist.tolist(), s=c1["s_D_mid_1.5"].tolist(),
                       GM1=float(c1["GM1_1.5"].mean()), FS=float(c1.rz_15.mean()))
    # radial profile of predictions for the median counter-rotator
    vc, re = float(t.Vc_star.median()), float(t.Re_kpc.median())
    radii = np.linspace(0.5, 4.0, 15)
    prof = {k: [] for k in MODELS}
    for m in radii:
        prof["FS"].append(pred_rz(fs, vc, m * re))
        prof["GM1"].append(predict(gm, vc, m * re, "n1"))
        prof["GM2"].append(predict(gm, vc, m * re, "n2"))
    out["profile"] = dict(Vc=vc, Re_kpc=re, radii=radii.tolist(), **{k: list(map(float, v)) for k, v in prof.items()})
    (ROOT / "results" / "paper_numbers.json").write_text(json.dumps(out, indent=2, default=float))
    t.to_csv(ROOT / "results" / "paper_manga_cr.csv")
    figures(t, out, c1)


def figures(t, out, c1):
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})
    col = {"FS": "#1b9e77", "GM1": "#d95f02", "GM2": "#7570b3"}
    # Fig. 1: predicted asymmetry against radius
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    p = out["profile"]
    for k, lab in MODELS.items():
        ax.plot(p["radii"], p[k], color=col[k], lw=1.4, label=lab)
    for tag, x in (("1.0", 1.0), ("1.5", 1.5)):
        a = out[f"R{tag}"]["intercepts"]["GM1_D_mid"]
        ax.errorbar(x, a["a"], yerr=a["sd_a"], fmt="o", color="k", ms=4, capsize=2,
                    label="MaNGA zero-twist intercept" if tag == "1.0" else None)
    ax.axhline(0, color="0.5", lw=0.6, ls=":")
    ax.set_xlabel(r"$R/R_{\rm e}$")
    ax.set_ylabel(r"$\Delta V = v_{\rm pro}-|v_{\rm retro}|$ (km s$^{-1}$)")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    ax.set_xlim(0.4, 4.1)
    fig.tight_layout()
    fig.savefig(FIG / "fig_predictions.pdf")
    plt.close(fig)
    # Fig. 2: slowdown against twist
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
    for ax, tag in zip(axes, ("1.0", "1.5")):
        m = np.isfinite(t[f"s_D_mid_{tag}"]) & np.isfinite(t.twist)
        x, y = t.twist[m].values.astype(float), t[f"s_D_mid_{tag}"][m].values
        ax.scatter(x, y, s=9, c="0.35", lw=0, label="MaNGA (exploratory)")
        a, b = theil_sen(x, y)
        xx = np.linspace(0, 180, 50)
        ax.plot(xx, a[0] + b[0] * xx, color="k", lw=1)
        for k in MODELS:
            ax.axhline(out[f"R{tag}"]["pred"][k], color=col[k], lw=1.1, ls="--",
                       label=MODELS[k] + " (prediction)")
        if tag == "1.5":
            ax.scatter(out["SAMI"]["twist"], out["SAMI"]["s"], s=22, marker="D",
                       facecolor="none", edgecolor="#e7298a", lw=0.9, label="SAMI (registered, N=5)")
        ax.axhline(0, color="0.5", lw=0.6, ls=":")
        ax.set_xlabel("gas kinematic twist (deg)")
        ax.set_title(rf"$R = {tag}\,R_{{\rm e}}$", fontsize=8)
        ax.set_xlim(-5, 185)
        ax.set_ylim(-230, 270)
    axes[0].set_ylabel(r"slowdown $s$ (km s$^{-1}$)")
    axes[1].legend(frameon=False, fontsize=6, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "fig_twist.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
