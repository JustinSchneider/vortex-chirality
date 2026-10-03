#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Co- and counter-rotating circular geodesics in the "strong gravitomagnetic"
low-energy GR limit of Astesiano, Ruggiero & Re (2026, arXiv:2606.18868).

Their metric (Eq. 1) has a Newtonian-like potential Phi and a dragging
function psi(r, z). Circular geodesics at constant (r, z) obey their Eq. (4),
written in the convention where Phi has the sign opposite to the Newtonian
potential (their footnote 1), so that  d_r Phi = -V0^2 / r  in the plane, with
V0 the Newtonian (baryonic) circular speed:

    (V/r) d_r psi = d_r Phi + V^2/r + (psi/r) d_r(psi/r).

Because Eq. (4) is a geodesic equation, it holds for any freely orbiting test
body with its own V, including gas counter-rotating with respect to the
stars. The script:
  1. solves Eq. (4) for V and recovers their Eq. (17) for psi = C0 r;
  2. derives the general branch asymmetry  V+ - |V-| = d_r psi;
  3. expresses C0 for the solution psi = C0 r in terms of observables:
        C0 = (V_obs^2 - V0^2) / V_obs,
     i.e. the predicted co/counter speed difference equals the mass
     discrepancy in velocity units (independent of which branch the stars
     occupy);
  4. checks the limits C0 -> 0 (Newton) and V0 -> 0 (V+ -> C0, their Eq. 18);
  5. generalises to psi = K r^n: the asymmetry is n K r^(n-1), and for
     1 <= n <= 2 it is at least the velocity discrepancy V_obs - V0, so
     the prediction does not hinge on the specific choice psi = C0 r.

Exit code 0 = all assertions passed.
"""
import sys

import sympy as sp

r, V, V0, C0, Vobs = sp.symbols("r V V0 C0 V_obs", positive=True)
psi = sp.Function("psi")(r)

ok = True


def check(name, cond):
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok &= bool(cond)


# Eq. (4) in the plane, with d_r Phi = -V0^2/r
eq4 = sp.Eq((V / r) * sp.diff(psi, r),
            -V0**2 / r + V**2 / r + (psi / r) * sp.diff(psi / r, r))

# 1. psi = C0 r  ->  Eq. (17)
e = eq4.subs(psi, C0 * r).doit()
roots = sp.solve(sp.simplify(e.lhs - e.rhs), V)
vp = (C0 + sp.sqrt(C0**2 + 4 * V0**2)) / 2
vm = (C0 - sp.sqrt(C0**2 + 4 * V0**2)) / 2
check("psi = C0 r reproduces Eq. (17), both branches",
      {sp.simplify(x) for x in roots} == {sp.simplify(vp), sp.simplify(vm)})
check("V- < 0: the second branch is counter-rotating",
      sp.simplify(vm.subs({C0: 1, V0: 1})) < 0)

# 2. general psi: quadratic in V; the sum of the roots is d_r psi
quad = sp.expand(r * (eq4.lhs - eq4.rhs))
a2, a1 = quad.coeff(V, 2), quad.coeff(V, 1)
sum_roots = sp.simplify(-a1 / a2)
check("V+ + V- = d_r psi for any psi, so V+ - |V-| = d_r psi",
      sp.simplify(sum_roots - sp.diff(psi, r)) == 0)
check("for psi = C0 r: V+ - |V-| = C0",
      sp.simplify(vp - (-vm) - C0) == 0)

# 3. C0 from the observed co-rotating speed and the baryonic speed
c0_obs = sp.solve(sp.Eq(vp.subs(C0, sp.Symbol("c")), Vobs), sp.Symbol("c"))
check("C0 = (V_obs^2 - V0^2) / V_obs",
      len(c0_obs) == 1 and sp.simplify(c0_obs[0] - (Vobs**2 - V0**2) / Vobs) == 0)

# 4. limits
check("C0 -> 0: V+ = V0 (Newtonian)", sp.simplify(vp.subs(C0, 0) - V0) == 0)
check("V0 -> 0: V+ -> C0 (flat curve, their Eq. 18)", sp.limit(vp, V0, 0) == C0)
check("counter-rotating speed |V-| -> 0 as V0 -> 0", sp.limit(-vm, V0, 0) == 0)

# 5. general power law psi = K r^n (in the plane), with x = K r^(n-1):
#    the asymmetry is n x, and for 1 <= n <= 2 it is never smaller than the
#    plain velocity discrepancy V_obs - V0 (n = 2 is rigid dragging, which
#    reproduces flowing-space Formulation A: V+ = V0 + x, V- = x - V0).
n, K, x = sp.symbols("n K x", positive=True)
en = eq4.subs(psi, K * r**n).doit()
qn = sp.expand(sp.simplify(r * (en.lhs - en.rhs)))
qn = sp.simplify(qn.subs(K, x / r**(n - 1)))
b2, b1 = sp.Poly(qn, V).all_coeffs()[0], sp.Poly(qn, V).all_coeffs()[1]
check("psi = K r^n: V+ - |V-| = n K r^(n-1)", sp.simplify(-b1 / b2 - n * x) == 0)
vpn = sp.solve(qn, V)
vpn = [v for v in vpn if sp.simplify(v.subs({n: 1.5, x: 1, V0: 1})) > 0][0]
worst = min(float((n * x - (vpn - V0)).subs({n: nn, x: xx, V0: 1}))
            for nn in (1, 1.25, 1.5, 1.75, 2) for xx in (0.05, 0.2, 0.5, 1, 2, 5))
check("1 <= n <= 2: V+ - |V-| >= V_obs - V0 (grid check)", worst >= -1e-12)
check("n = 2: V+ = V0 + K r (flowing-space Formulation A)",
      sp.simplify(vpn.subs(n, 2) - (V0 + x)) == 0)

# Numerical illustration: V_obs = 200, V0 = 150 km/s
num = ((Vobs**2 - V0**2) / Vobs).subs({Vobs: 200, V0: 150})
print(f"  example: V_obs = 200, V_bar = 150 km/s  ->  predicted |V+| - |V-| = {float(num):.1f} km/s")

sys.exit(0 if ok else 1)
