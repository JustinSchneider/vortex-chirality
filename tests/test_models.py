import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from src import models

ROOT = Path(__file__).resolve().parents[1]
R = np.linspace(0.1, 30, 300)
VB = 80 * np.sqrt(R / (R + 2.0))      # smooth rising baryonic curve


def test_no_medium_recovers_newton():
    v = models.v_flowing_space(R, VB, u=0.0, shear=0.0)
    np.testing.assert_allclose(v, VB)
    np.testing.assert_allclose(models.v_force_field(R, VB, 0.0), VB)


def test_solid_body_flow_is_exactly_additive():
    Om = 7.0
    v = models.v_flowing_space(R, VB, u=Om * R, shear=0.0)
    np.testing.assert_allclose(v, VB + Om * R, rtol=1e-12)


def test_force_field_weak_limit_is_fc25_linear():
    Om = 1e-3
    v = models.v_force_field(R, VB, Om)
    np.testing.assert_allclose(v, VB + Om * R, rtol=1e-6)


def test_circular_orbit_condition_A():
    """v from the model satisfies w^2 + S R w - Vb^2 = 0."""
    om, Rt = 25.0, 3.0
    for retro in (False, True):
        v = models.v_flowing_space_rt(R, VB, om, Rt, retrograde=retro)
        u = models.rt_flow(R, om, Rt)
        S = models.rt_flow_shear(R, om, Rt)
        w = v - u
        np.testing.assert_allclose(w**2 + S * R * w - VB**2, 0, atol=1e-8)


def test_shear_and_vorticity_match_finite_differences():
    om, Rt = 25.0, 3.0
    r = np.linspace(0.1, 30, 30000)
    u = models.rt_flow(r, om, Rt)
    S_num = u / r - np.gradient(u, r)
    zeta_num = np.gradient(r * u, r) / r
    sl = slice(5, -5)
    np.testing.assert_allclose(models.rt_flow_shear(r, om, Rt)[sl],
                               S_num[sl], rtol=1e-5)
    np.testing.assert_allclose(models.rt_flow_vorticity(r, om, Rt)[sl],
                               zeta_num[sl], rtol=1e-5)


@pytest.mark.parametrize("om,Rt", [(25.0, 3.0), (10.0, 10.0), (60.0, 1.0)])
def test_asymmetry_equals_R_times_vorticity(om, Rt):
    vp = models.v_flowing_space_rt(R, VB, om, Rt)
    vr = models.v_flowing_space_rt(R, VB, om, Rt, retrograde=True)
    # Delta V = v_pro - |v_ret|, valid where the retrograde orbit is retrograde
    mask = vr < 0
    np.testing.assert_allclose((vp + vr)[mask],
                               (R * models.rt_flow_vorticity(R, om, Rt))[mask],
                               rtol=1e-10)


def test_force_field_asymmetry_is_2_omega_R():
    Om = 9.0
    dv = models.v_force_field(R, VB, Om) - models.v_force_field(
        R, VB, Om, retrograde=True)
    np.testing.assert_allclose(dv, 2 * Om * R, rtol=1e-12)


def test_retrograde_branch_is_prograde_with_reversed_flow():
    """A retrograde orbit in flow +u mirrors a prograde orbit in flow -u."""
    om, Rt = 25.0, 3.0
    vr = models.v_flowing_space_rt(R, VB, om, Rt, retrograde=True)
    vp_rev = models.v_flowing_space_rt(R, VB, -om, Rt)
    np.testing.assert_allclose(-vr, vp_rev, rtol=1e-12)


def test_chirality_estimator_recovers_injected_flow():
    """Build two drift-lagged counter-rotating disks in a solid-body flow."""
    Vb, u, k = 150.0, 20.0, 2.0
    s1, s2 = 30.0, 70.0
    # medium-frame speeds w satisfy w^2 + k s^2 = Vb^2
    w1, w2 = np.sqrt(Vb**2 - k * s1**2), np.sqrt(Vb**2 - k * s2**2)
    V1, V2 = w1 + u, w2 - u          # observed speeds (both positive)
    est = models.chirality_estimator(V1, V2, s1, s2, k)
    assert est == pytest.approx(2 * u, rel=1e-12)


def test_chirality_estimator_mixed_stars_and_cold_gas():
    """Stars (drift-lagged) and cold gas (k = 0) counter-rotating."""
    Vb, u, k = 150.0, 15.0, 2.5
    s_star = 60.0
    w_star = np.sqrt(Vb**2 - k * s_star**2)
    V_star, V_gas = w_star + u, Vb - u
    est = models.chirality_estimator(V_star, V_gas, s_star, 0.0, k, k2=0.0)
    assert est == pytest.approx(2 * u, rel=1e-12)


@pytest.mark.parametrize("script", sorted((ROOT / "derivations").glob("0*.py")))
def test_derivation_scripts_pass(script):
    r = subprocess.run([sys.executable, str(script)], capture_output=True,
                       text=True, encoding="utf-8", cwd=ROOT)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
