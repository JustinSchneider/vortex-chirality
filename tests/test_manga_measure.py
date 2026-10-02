"""Synthetic validation of the MaNGA stage-2 pipeline (PREREG f.9).

A mock galaxy has baryonic circular speed Vb(R), stars and counter-rotating
gas with known dispersion profiles, embedded in a solid-body flow u = Omega R.
Formulation A then fixes their observed speeds exactly (see
derivations/01_force_laws.py section C), and the pipeline must recover
D = -R*zeta = -2*Omega*R at 1 R_e to within 10%.
"""

import numpy as np
import pytest

from src import manga_measure as mm

RNG = np.random.default_rng(7)
RE = 6.0            # arcsec
INC, PA_S = 55.0, 30.0
PA_G = PA_S + 180.0


def mock(Omega, sig0_s=90.0, sig0_g=25.0, noise=4.0, n=64, vr_gas=0.0):
    yy, xx = np.mgrid[0:n, 0:n]
    x, y = (xx - (n - 1) / 2) * 0.5, (yy - (n - 1) / 2) * 0.5
    h_R = RE / mm.EXP_RE_TO_H
    h_s2_s, h_s2_g = 2.0 * h_R, 3.0 * h_R

    def comp(pa, sig0, h_s2, sense, vr=0.0):
        R, phi = mm.disk_coords(x, y, pa, INC)
        Vb = 220.0 * R / np.sqrt(R**2 + (0.4 * RE) ** 2)
        sig = sig0 * np.exp(-R / (2 * h_s2))            # sigma^2 e-folds at h_s2
        k = np.clip(R / h_R + R / h_s2 - 0.5, *mm.K_CLIP)
        w = np.sqrt(np.clip(Vb**2 - k * sig**2, 0, None))
        V = w + sense * Omega * R                       # observed speed (>0)
        vlos = (V * np.cos(phi) + vr * np.sin(phi)) * np.sin(np.radians(INC))
        vlos += RNG.normal(0, noise, vlos.shape)
        return vlos, np.full_like(vlos, noise), sig

    vs, evs, ss = comp(PA_S, sig0_s, h_s2_s, +1)   # stars co-rotate with flow
    vg, evg, sg = comp(PA_G, sig0_g, h_s2_g, -1, vr_gas)   # gas counter-rotates
    maps = dict(vs=vs, evs=evs, ss=ss, vg=vg, evg=evg, sg=sg)
    return maps, x, y


@pytest.mark.parametrize("Omega_kms_arcsec", [0.0, 3.0, 6.0])
def test_recovers_injected_asymmetry(Omega_kms_arcsec):
    maps, x, y = mock(Omega_kms_arcsec)
    rows = mm.measure(maps, x, y, PA_S, PA_G, INC, RE)
    r1 = next(r for r in rows if r["r_re"] == 1.0)
    truth = -2 * Omega_kms_arcsec * RE
    tol = max(0.10 * abs(truth), 3.0)
    assert r1["D_mid"] == pytest.approx(truth, abs=tol)


def test_null_flow_gives_D_near_zero_at_all_rings():
    maps, x, y = mock(0.0)
    rows = mm.measure(maps, x, y, PA_S, PA_G, INC, RE)
    for r in rows:
        if np.isfinite(r["D_mid"]):
            assert abs(r["D_mid"]) < 8.0, r


def test_disk_coords_receding_axis():
    R, phi = mm.disk_coords(np.array([np.sin(np.radians(30)) * 2]),
                            np.array([np.cos(np.radians(30)) * 2]), 30.0, 60.0)
    assert R[0] == pytest.approx(2.0)
    assert phi[0] == pytest.approx(0.0, abs=1e-12)


def test_radial_inflow_recovered_and_does_not_bias_rotation():
    m0, x, y = mock(3.0)
    m1, _, _ = mock(3.0, vr_gas=-25.0)
    r0 = next(r for r in mm.measure(m0, x, y, PA_S, PA_G, INC, RE) if r["r_re"] == 1.0)
    r1 = next(r for r in mm.measure(m1, x, y, PA_S, PA_G, INC, RE) if r["r_re"] == 1.0)
    assert r1["VRg"] == pytest.approx(25.0, abs=3.0)
    assert r0["VRg"] < 3.0
    assert r1["D_mid"] == pytest.approx(r0["D_mid"], abs=3.0)
