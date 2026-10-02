#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Dimensional analysis (Buckingham Pi) for the scaling laws the vortical-
vacuum and seed-first hypotheses are allowed to take.

Every quantity is represented by its exponent vector over (M, L, T). A
product of quantities is dimensionless iff the exponent-weighted sum of
their vectors is zero, so the dimensionless groups span the nullspace of
the dimension matrix. This lets the script *derive* the Pi groups rather
than assert them.

Sections:
  1. Pi groups for the RT parameters (R_t, R_d, M_bar, M_BH, G, a0).
  2. What a theory with NO disk length scale predicts (pure-a0 scaling),
     and how existing data reject it.
  3. Units of the coupling kappa in curl(Omega) = kappa * rho * v, and the
     natural constants with those units.
  4. The unique combination that gives the observed omega magnitude.
  5. Dimensional consistency of the force laws from 01_force_laws.py.

Exit code 0 = all assertions passed.
"""

import sys

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Dimension exponent vectors over (M, L, T)
DIM = {
    "R_t":   (0, 1, 0),
    "R_d":   (0, 1, 0),
    "R":     (0, 1, 0),
    "M_bar": (1, 0, 0),
    "M_BH":  (1, 0, 0),
    "G":     (-1, 3, -2),
    "a0":    (0, 1, -2),
    "c":     (0, 1, -1),
    "V":     (0, 1, -1),
    "omega": (0, 0, -1),
    "rho":   (1, -3, 0),
    "H0":    (0, 0, -1),
}


def dim_of(**powers):
    """Exponent vector of prod(name**power)."""
    vec = sp.zeros(3, 1)
    for name, p in powers.items():
        vec += sp.Rational(p) * sp.Matrix(DIM[name])
    return vec


def is_dimensionless(**powers):
    return dim_of(**powers) == sp.zeros(3, 1)


def check(label, cond):
    print(f"    [{'PASS' if cond else 'FAIL'}] {label}")
    if not cond:
        sys.exit(1)


def pi_groups(names):
    """Return a basis of dimensionless groups for the given variables."""
    M = sp.Matrix([list(DIM[n]) for n in names]).T  # 3 x n
    basis = M.nullspace()
    groups = []
    for vec in basis:
        # Clear denominators so exponents are integers.
        lcm = sp.ilcm(*[sp.fraction(x)[1] for x in vec])
        vec = vec * lcm
        groups.append({n: int(e) for n, e in zip(names, vec) if e != 0})
    return M.rank(), groups


def fmt(group):
    return " * ".join(f"{n}^{e}" if e != 1 else n for n, e in group.items())


print("=" * 72)
print("1. Buckingham Pi for the RT parameters")
print("=" * 72)
names = ["R_t", "R_d", "M_bar", "M_BH", "G", "a0"]
rank, groups = pi_groups(names)
print(f"\n    {len(names)} variables, dimension-matrix rank {rank}"
      f" -> {len(names) - rank} independent groups. A nullspace basis:")
for g in groups:
    print(f"      {fmt(g)}")
check("6 variables, rank 3 -> exactly 3 Pi groups",
      len(names) - rank == 3 and len(groups) == 3)

# The physically readable basis used in docs/theory.md:
Pi1 = dict(R_t=1, R_d=-1)
Pi2 = dict(G=1, M_bar=1, a0=-1, R_d=-2)
Pi3 = dict(M_BH=1, M_bar=-1)
for label, g in [("Pi1 = R_t/R_d", Pi1),
                 ("Pi2 = G M_bar / (a0 R_d^2)  [= Sigma_d / Sigma_dagger]", Pi2),
                 ("Pi3 = M_BH / M_bar", Pi3)]:
    check(f"{label} is dimensionless", is_dimensionless(**g))

# Independence of the readable basis.
mat = sp.Matrix([[g.get(n, 0) for n in names] for g in (Pi1, Pi2, Pi3)])
check("Pi1, Pi2, Pi3 are independent (rank 3)", mat.rank() == 3)
print("""
    => The most general allowed law is
         R_t / R_d            = F(Pi2, Pi3)
         V_sat^4/(G M_bar a0) = H(Pi2, Pi3)
       H_seed predicts dF/dPi3 != 0 at fixed Pi2.""")

# Velocity group used for the BTFR-like relation.
check("V_sat^4 / (G M_bar a0) is dimensionless",
      is_dimensionless(V=4, G=-1, M_bar=-1, a0=-1))


print("\n" + "=" * 72)
print("2. A theory with no disk length scale (only G, M_bar, a0)")
print("=" * 72)
G, M, a0 = sp.symbols("G M a_0", positive=True)
# With only G, M, a0 the unique length, speed and rate are:
L_M = sp.sqrt(G * M / a0)          # MOND radius
V_M = (G * M * a0) ** sp.Rational(1, 4)
W_M = sp.simplify(V_M / L_M)
check("sqrt(G M / a0) has dimensions of length",
      dim_of(G=sp.Rational(1, 2), M_bar=sp.Rational(1, 2),
             a0=-sp.Rational(1, 2)) == sp.Matrix(DIM["R"]))
check("(G M a0)^(1/4) has dimensions of speed",
      dim_of(G=sp.Rational(1, 4), M_bar=sp.Rational(1, 4),
             a0=sp.Rational(1, 4)) == sp.Matrix(DIM["V"]))
print(f"\n    R_t   ~ {L_M}           -> R_t   ∝ M^{sp.Rational(1, 2)}")
print(f"    V_sat ~ {V_M}  -> V_sat ∝ M^{sp.Rational(1, 4)}")
print(f"    omega ~ {W_M}")
check("pure-a0 theory forces omega * V_sat = a0 (a constant)",
      sp.simplify(W_M * V_M - a0) == 0)
print("""
    In Pi-group form this is R_t/R_d ∝ Pi2^(1/2): a pure-a0 theory still
    makes F depend on Pi2, with a fixed exponent of 1/2. A pure disk-scale
    theory (R_t ∝ R_d) gives exponent 0. The exponent is measured, not
    assumed, in analysis/s2_scaling_law.py (tests S2a, S2b).""")


print("\n" + "=" * 72)
print("3. Coupling for a current-sourced vorticity: curl(Omega) = kappa rho v")
print("=" * 72)
# [curl Omega] = [Omega]/L ; [rho v] = M L^-2 T^-1 ; kappa = ratio.
kappa_dim = (dim_of(omega=1, R=-1) - dim_of(rho=1, V=1))
print(f"\n    [kappa] exponents (M, L, T) = {tuple(kappa_dim)}  -> L / M")
check("[kappa] = L/M", kappa_dim == sp.Matrix([-1, 1, 0]))
check("G/c^2 has dimensions L/M (the GR gravitomagnetic coupling)",
      dim_of(G=1, c=-2) == kappa_dim)
check("G/V^2 for ANY speed V also has dimensions L/M",
      dim_of(G=1, V=-2) == kappa_dim)
print("""
    The natural form is G/c_v^2 for some characteristic speed c_v.
    GR has c_v = c; matching the observed omega needs c_v << c (see
    03_magnitude_gates.py for the number).""")


print("\n" + "=" * 72)
print("4. The combination that sets the omega magnitude")
print("=" * 72)
# Find exponents (x, y) with a0^x V^y having dimensions of 1/T.
x, y = sp.symbols("x y")
eqs = list(x * sp.Matrix(DIM["a0"]) + y * sp.Matrix(DIM["V"])
           - sp.Matrix(DIM["omega"]))
sol = sp.solve(eqs, [x, y], dict=True)
print(f"\n    a0^x V^y with dimensions 1/T: {sol}")
check("unique solution Omega ~ a0 / V", sol == [{x: 1, y: -1}])
check("c H0 has dimensions of acceleration (a0 ≈ c H0 / 2pi)",
      dim_of(c=1, H0=1) == sp.Matrix(DIM["a0"]))


print("\n" + "=" * 72)
print("5. Dimensional consistency of the force laws")
print("=" * 72)
# v^2/R, Vb^2/R, Omega v, and S R w all must be accelerations.
acc = sp.Matrix(DIM["a0"])
check("v^2 / R is an acceleration", dim_of(V=2, R=-1) == acc)
check("Omega * v is an acceleration", dim_of(omega=1, V=1) == acc)
check("Omega * R is a speed", dim_of(omega=1, R=1) == sp.Matrix(DIM["V"]))
check("R * zeta (Delta V) is a speed",
      dim_of(R=1, omega=1) == sp.Matrix(DIM["V"]))
check("omega R / (1 + R/R_t) is a speed (R/R_t dimensionless)",
      is_dimensionless(R=1, R_t=-1)
      and dim_of(omega=1, R=1) == sp.Matrix(DIM["V"]))

print("\nAll assertions passed.")
