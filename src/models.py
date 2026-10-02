"""Rotation-curve models derived in derivations/01_force_laws.py.

Units: R in kpc, speeds in km/s, omega / Omega in km/s/kpc.

Formulation A ("flowing space"): inertia is relative to a medium moving at
u(R) azimuthally. A circular orbit satisfies w^2 + S R w - Vb^2 = 0 with
w = v - u and shear S = u/R - du/dR, so
    v = u - S R/2 +/- sqrt(S^2 R^2/4 + Vb^2).

Formulation B ("force field"): a = g + 2 v x Omega, so
    v = +/- Omega R + sqrt(Omega^2 R^2 + Vb^2).

Both predict a prograde/retrograde asymmetry Delta V = R * zeta(R).
"""

import numpy as np


def rt_flow(R, omega, Rt):
    """RT profile read as a flow speed: u = omega R / (1 + R/Rt)."""
    R = np.asarray(R, dtype=float)
    return omega * R / (1.0 + R / Rt)


def rt_flow_shear(R, omega, Rt):
    """Shear S = u/R - du/dR for the RT flow (km/s/kpc).

    u/R = omega/(1+x), du/dR = omega/(1+x)^2 with x = R/Rt, so
    S = omega x / (1+x)^2.
    """
    x = np.asarray(R, dtype=float) / Rt
    return omega * x / (1.0 + x) ** 2


def rt_flow_vorticity(R, omega, Rt):
    """Vorticity zeta = (1/R) d(R u)/dR of the RT flow (km/s/kpc)."""
    R = np.asarray(R, dtype=float)
    return omega * Rt * (R + 2 * Rt) / (R + Rt) ** 2


def v_flowing_space(R, v_bary, u, shear, retrograde=False):
    """Observed speed in Formulation A (signed; retrograde orbits < 0)."""
    R = np.asarray(R, dtype=float)
    root = np.sqrt((shear * R) ** 2 / 4.0 + np.asarray(v_bary) ** 2)
    sign = -1.0 if retrograde else 1.0
    return u - shear * R / 2.0 + sign * root


def v_flowing_space_rt(R, v_bary, omega, Rt, retrograde=False):
    """Formulation A with the RT flow profile (2 parameters: omega, Rt)."""
    return v_flowing_space(R, v_bary, rt_flow(R, omega, Rt),
                           rt_flow_shear(R, omega, Rt), retrograde)


def v_force_field(R, v_bary, Omega, retrograde=False):
    """Observed speed magnitude in Formulation B."""
    R = np.asarray(R, dtype=float)
    sign = -1.0 if retrograde else 1.0
    return sign * Omega * R + np.sqrt((Omega * R) ** 2
                                      + np.asarray(v_bary) ** 2)


def v_rt_additive(R, v_bary, omega, Rt):
    """Original Rational Taper: V = Vb + omega R / (1 + R/Rt)."""
    return np.asarray(v_bary) + rt_flow(R, omega, Rt)


def chirality_estimator(V1, V2, s1, s2, k, k2=None):
    """R*zeta from two counter-rotating disks (solid-body flow).

    V1, V2: rotation speeds (both positive) of the two disks at the same
    radius; s1, s2: their velocity dispersions; k: asymmetric-drift factor
    (v_c^2 - v^2 = k sigma^2) for disk 1, and for disk 2 unless k2 is given
    (use k = 0 for cold gas). Positive result means the flow co-rotates
    with disk 1. Derived in derivations/01_force_laws.py, section C.
    """
    k2 = k if k2 is None else k2
    return (V1 ** 2 - V2 ** 2 + k * s1 ** 2 - k2 * s2 ** 2) / (V1 + V2)
