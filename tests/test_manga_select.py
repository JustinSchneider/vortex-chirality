"""The selection stage must recover a known receding-side PA, in the same
convention that src/manga_measure.disk_coords uses (degrees from +y towards
+x). This is the test that would have caught the pafit convention bug."""

import numpy as np
import pytest

from analysis.g1_manga_select import grid_xy, receding_pa
from src import manga_measure as mm


def velfield(pa_true, inc=55.0, vmax=150.0, n=40):
    x, y = grid_xy((n, n))
    R, phi = mm.disk_coords(x, y, pa_true, inc)
    v = vmax * np.tanh(R / 3.0) * np.cos(phi) * np.sin(np.radians(inc))
    return x.ravel(), y.ravel(), v.ravel()


@pytest.mark.parametrize("pa_true", [10, 60, 100, 170, 200, 250, 300, 345])
def test_receding_pa_recovered(pa_true):
    x, y, v = velfield(pa_true)
    pa, err = receding_pa(x, y, v, np.full_like(v, 5.0))
    d = abs(pa - pa_true) % 360
    assert min(d, 360 - d) < 3.0, (pa_true, pa)


@pytest.mark.parametrize("pa_s", [20, 75, 130])
def test_counter_rotation_gives_dpa_180(pa_s):
    xs, ys, vs = velfield(pa_s)
    xg, yg, vg = velfield(pa_s + 180)
    ps, _ = receding_pa(xs, ys, vs, np.full_like(vs, 5.0))
    pg, _ = receding_pa(xg, yg, vg, np.full_like(vg, 5.0))
    d = abs(pg - ps) % 360
    assert min(d, 360 - d) > 175
