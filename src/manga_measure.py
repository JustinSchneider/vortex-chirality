"""MaNGA G1 stage 2: ring rotation curves, dispersions, drift factors and D.

Implements PREREG amendment (f). All functions take plain arrays so they can
be validated on synthetic maps (tests/test_manga_measure.py) before any real
galaxy is measured.

Conventions: x, y are sky offsets (arcsec) from the map centre; position
angles are receding-side PAs in degrees measured from +y towards +x (same
convention as analysis/g1_manga_select.receding_pa).
"""

import numpy as np

RING_CENTRES = (0.5, 0.75, 1.0, 1.25, 1.5)   # in units of R_e
RING_WIDTH = 0.25
MIN_RING_SPAX = 20
K_CLIP = (0.0, 4.0)
EXP_RE_TO_H = 1.678                           # R_e = 1.678 h for an exponential


def disk_coords(x, y, pa_deg, inc_deg):
    """In-plane radius (arcsec) and azimuth (rad, 0 on the receding axis)."""
    th = np.radians(pa_deg)
    # major-axis coordinate (towards receding PA) and minor-axis coordinate
    xm = x * np.sin(th) + y * np.cos(th)
    ym = -x * np.cos(th) + y * np.sin(th)
    yd = ym / np.cos(np.radians(inc_deg))
    return np.hypot(xm, yd), np.arctan2(yd, xm)


def ring_fit(v, verr, R, phi, r_lo, r_hi, inc_deg):
    """Harmonic fit v = Vsys + c1 cos(phi) + s1 sin(phi) in one ring.

    Returns (V_rot, err, n, resid_rms) or None if the ring is invalid.
    """
    m = (R >= r_lo) & (R < r_hi) & np.isfinite(v)
    n = int(m.sum())
    if n < MIN_RING_SPAX:
        return None
    quad = np.floor(((phi[m] + np.pi) % (2 * np.pi)) / (np.pi / 2)).astype(int)
    if len(np.unique(quad)) < 3:
        return None
    A = np.column_stack([np.ones(n), np.cos(phi[m]), np.sin(phi[m])])
    w = 1.0 / np.clip(verr[m], 1.0, None)
    coef, *_ = np.linalg.lstsq(A * w[:, None], v[m] * w, rcond=None)
    resid = v[m] - A @ coef
    cov = np.linalg.pinv((A * w[:, None]).T @ (A * w[:, None]))
    chi2r = np.sum((resid * w) ** 2) / max(n - 3, 1)
    sin_i = np.sin(np.radians(inc_deg))
    V = coef[1] / sin_i
    eV = np.sqrt(cov[1, 1] * max(chi2r, 1.0)) / sin_i
    # s1 term = projected radial (in/outflow) speed; its sign depends on the
    # unknown near side, so only |V_R| is meaningful.
    VR = coef[2] / sin_i
    return float(V), float(eV), n, float(np.sqrt(np.mean(resid ** 2))), float(VR)


def ring_median(q, R, r_lo, r_hi):
    m = (R >= r_lo) & (R < r_hi) & np.isfinite(q)
    return float(np.median(q[m])) if m.sum() >= MIN_RING_SPAX else np.nan


def rotation_profile(v, verr, sig, x, y, pa_deg, inc_deg, Re):
    """V_rot and sigma at each ring centre (R in units of R_e)."""
    R, phi = disk_coords(x, y, pa_deg, inc_deg)
    out = []
    for c in RING_CENTRES:
        lo, hi = (c - RING_WIDTH / 2) * Re, (c + RING_WIDTH / 2) * Re
        f = ring_fit(v, verr, R, phi, lo, hi, inc_deg)
        out.append(dict(r_re=c,
                        V=f[0] if f else np.nan, eV=f[1] if f else np.nan,
                        n=f[2] if f else 0,
                        resid=f[3] if f else np.nan,
                        VR=abs(f[4]) if f else np.nan,
                        sigma=ring_median(sig, R, lo, hi)))
    return out


def drift_k(profile, Re):
    """k(R) = R/h_R + R/h_{sigma^2} - 1/2 at each ring, clipped to [0, 4]."""
    r = np.array([p["r_re"] for p in profile]) * Re
    s2 = np.array([p["sigma"] for p in profile]) ** 2
    good = np.isfinite(s2) & (s2 > 0)
    h_R = Re / EXP_RE_TO_H
    if good.sum() >= 3:
        slope = np.polyfit(r[good], np.log(s2[good]), 1)[0]   # = -1/h_s2
        inv_h_s2 = -slope
    else:
        inv_h_s2 = 0.0
    k = r / h_R + r * inv_h_s2 - 0.5
    return np.clip(k, *K_CLIP)


def d_statistic(Vg, sg, kg, Vs, ss, ks):
    """D = gas minus drift-corrected stellar circular speed (km/s)."""
    return (Vg ** 2 + kg * sg ** 2 - Vs ** 2 - ks * ss ** 2) / (Vg + Vs)


def measure(maps, x, y, pa_star, pa_gas, inc, Re_arcsec):
    """Full stage-2 measurement for one galaxy from pre-masked maps.

    maps: dict with arrays 'vs','evs','ss' (stars) and 'vg','evg','sg' (gas);
    invalid spaxels must be NaN. Returns per-ring dict including D at
    k-1, k, k+1 (shared shift, PREREG f.4).
    """
    ps = rotation_profile(maps["vs"], maps["evs"], maps["ss"], x, y,
                          pa_star, inc, Re_arcsec)
    pg = rotation_profile(maps["vg"], maps["evg"], maps["sg"], x, y,
                          pa_gas, inc, Re_arcsec)
    ks, kg = drift_k(ps, Re_arcsec), drift_k(pg, Re_arcsec)
    rows = []
    for i, c in enumerate(RING_CENTRES):
        Vs, Vg = abs(ps[i]["V"]), abs(pg[i]["V"])
        row = dict(r_re=c, Vs=Vs, eVs=ps[i]["eV"], ss=ps[i]["sigma"],
                   Vg=Vg, eVg=pg[i]["eV"], sg=pg[i]["sigma"],
                   ks=ks[i], kg=kg[i], gas_resid=pg[i]["resid"],
                   VRg=pg[i]["VR"], VRs=ps[i]["VR"])
        for tag, dk in (("lo", -1.0), ("mid", 0.0), ("hi", 1.0)):
            row[f"D_{tag}"] = d_statistic(Vg, row["sg"], max(kg[i] + dk, 0),
                                          Vs, row["ss"], max(ks[i] + dk, 0))
        row["Vc_star"] = np.sqrt(Vs ** 2 + ks[i] * row["ss"] ** 2)
        rows.append(row)
    return rows
