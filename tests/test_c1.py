"""Checks for the confirmatory C1 pipeline on synthetic inputs."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis import c1_test  # noqa: E402
from analysis.c1_power import theil_sen  # noqa: E402
from src import sami  # noqa: E402


def _prim(a, b, n=40, sd=5.0, seed=1):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 90, n)
    s = a + b * t + rng.normal(0, sd, n)
    return pd.DataFrame({"twist": t, "s_D_mid_1.5": s})


def test_theil_sen_exact_line():
    x = np.arange(10.0)
    a, b = theil_sen(x, 3 + 0.5 * x)
    assert a[0] == pytest.approx(3) and b[0] == pytest.approx(0.5)


@pytest.mark.parametrize("a_true,expect", [(0.0, -2), (40.0, 2)])
def test_primary_recovers_intercept_and_decides(a_true, expect):
    r = c1_test.primary(_prim(a_true, 0.3), P=40.0, k="D_mid")
    assert r["a"] == pytest.approx(a_true, abs=4)
    assert r["b"] == pytest.approx(0.3, abs=0.1)
    assert r["code"] == expect
    assert r["f_ci95"][0] <= a_true / 40 <= r["f_ci95"][1]


def test_label_thresholds():
    assert c1_test.label(10)[0] == 2 and c1_test.label(5)[0] == 1
    assert c1_test.label(1)[0] == 0 and c1_test.label(0.2)[0] == -1
    assert c1_test.label(0.1)[0] == -2


def test_sami_find_file_unique(tmp_path):
    names = ["123_A_stellar-velocity_default_two-moment.fits",
             "123_A_stellar-velocity-dispersion_default_two-moment.fits",
             "123_A_gas-velocity_default_1-comp.fits",
             "123_A_gas-vdisp_default_1-comp.fits",
             "123_A_Halpha_default_1-comp.fits",
             "123_A_gas-velocity_default_recom-comp.fits",
             "1234_A_stellar-velocity_default_two-moment.fits"]
    files = [tmp_path / n for n in names]
    assert sami.find_file(123, "stellar_vel", files).name == names[0]
    assert sami.find_file(123, "stellar_sig", files).name == names[1]
    assert sami.find_file(123, "gas_vel", files).name == names[2]
    assert sami.find_file(123, "gas_sig", files).name == names[3]
    assert sami.find_file(123, "ha_flux", files).name == names[4]


def test_califa_overlap_galaxy_matches_manga():
    """Orientation check on a MaNGA-overlap galaxy (excluded from C1)."""
    from src import califa
    from analysis.g1_manga_select import receding_pa
    if not (califa.CACHE / "UGC08107_gas.npz").exists():
        pytest.skip("overlap galaxy not cached")
    r = califa.catalogue().loc["UGC08107"]
    maps, x, y, s_ok, g_ok = califa.galaxy_arrays("UGC08107", r.ra, r.dec, r.z)
    dec = np.zeros_like(s_ok)
    dec[::2, ::2] = True
    ps = receding_pa(x[s_ok & dec], y[s_ok & dec], maps["vs"][s_ok & dec], maps["evs"][s_ok & dec])[0]
    pg = receding_pa(x[g_ok & dec], y[g_ok & dec], maps["vg"][g_ok & dec], maps["evg"][g_ok & dec])[0]
    # MaNGA 11761-12705: both 309 deg
    for pa in (ps, pg):
        d = abs(pa - 309) % 360
        assert min(d, 360 - d) < 20


def test_sami_overlap_galaxy_matches_manga():
    """Orientation/format check on a MaNGA-overlap galaxy (excluded from C1)."""
    from analysis.g1_manga_select import receding_pa
    if not sami.available("517164"):
        pytest.skip("overlap galaxy not downloaded")
    maps, x, y, s_ok, g_ok = sami.galaxy_arrays("517164")
    dec = np.zeros_like(s_ok)
    dec[::2, ::2] = True
    ps = receding_pa(x[s_ok & dec], y[s_ok & dec], maps["vs"][s_ok & dec], maps["evs"][s_ok & dec])[0]
    pg = receding_pa(x[g_ok & dec], y[g_ok & dec], maps["vg"][g_ok & dec], maps["evg"][g_ok & dec])[0]
    # MaNGA 11754-1901: both 83 deg
    for pa in (ps, pg):
        d = abs(pa - 83) % 360
        assert min(d, 360 - d) < 20
