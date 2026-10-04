#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Symbolic derivation of rotation curves under the two formulations of the
vortical-vacuum hypothesis (H_vort), and of their chirality prediction.

Formulation A -- "flowing space" (river / frame-dragging picture):
    Inertia is defined relative to a medium that moves with velocity u(R)
    in the azimuthal direction:
        L_A = 1/2 |v - u|^2 - Phi_N
    (Unruh 1981 acoustic metric; Hamilton & Lisle 2008 river model.)

Formulation B -- "vortical force field" (gravitomagnetic, as in MGEM eq. 1):
    The medium exerts a velocity-dependent force with no centrifugal term:
        L_B = 1/2 v^2 - Phi_N + v . A,     (curl A)_z = B(R)

For each formulation the script:
  1. derives the radial Euler-Lagrange equation for a circular orbit,
  2. solves it for the observed (coordinate) speed of prograde and
     retrograde orbits,
  3. checks the limits (no medium -> Newton; weak medium -> FC25 linear),
  4. derives the prograde-minus-retrograde speed asymmetry
        Delta V(R) = v_pro - |v_retro|
     which is the decisive observable (Gate G1).

Exit code 0 = all assertions passed.
"""

import sys

import sympy as sp

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

t = sp.Symbol("t")
R_s, Vb, Om, w, v = sp.symbols("R V_b Omega w v", positive=True)
Rt, om = sp.symbols("R_t omega", positive=True)

Rf = sp.Function("R")(t)
phif = sp.Function("phi")(t)
Phi = sp.Function("Phi")
u = sp.Function("u")
A = sp.Function("A")


def radial_circular_equation(L):
    """Radial Euler-Lagrange equation evaluated on a circular orbit.

    Returns an expression that equals zero for a circular orbit with
    angular velocity phidot = v/R at radius R.
    """
    Rd, phid = sp.diff(Rf, t), sp.diff(phif, t)
    el = sp.diff(sp.diff(L, Rd), t) - sp.diff(L, Rf)
    el = el.subs({sp.diff(Rf, t, 2): 0})
    el = el.subs({phid: v / Rf, Rd: 0})
    el = el.subs(Rf, R_s).doit()
    # Newtonian gravity from baryons: dPhi/dR = V_b^2 / R
    el = el.subs(sp.Derivative(Phi(R_s), R_s), Vb**2 / R_s)
    return sp.simplify(el)


def check(label, expr):
    ok = sp.simplify(expr) == 0
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        print(f"           residual: {sp.simplify(expr)}")
        sys.exit(1)


print("=" * 72)
print("Formulation A: flowing space, L = 1/2 |v - u|^2 - Phi")
print("=" * 72)

Rd, phid = sp.diff(Rf, t), sp.diff(phif, t)
L_A = sp.Rational(1, 2) * (Rd**2 + (Rf * phid - u(Rf))**2) - Phi(Rf)
eq_A = radial_circular_equation(L_A)
print(f"\n[A1] Circular-orbit condition (=0):\n     {eq_A}")

# Write v = u + w, where w is the speed relative to the local medium.
U, Up = sp.symbols("u u_p")  # u(R) and du/dR at the orbit radius
eq_A_w = eq_A.subs(sp.Derivative(u(R_s), R_s), Up).subs(u(R_s), U)
eq_A_w = sp.expand(eq_A_w.subs(v, U + w) * R_s)
print(f"[A2] In terms of w = v - u (times R):\n     {sp.factor(eq_A_w)} = 0")

# Shear term S = u/R - du/dR. The condition is w^2 + S R w - Vb^2 = 0.
S = sp.Symbol("S")
target = w**2 + (U / R_s - Up) * R_s * w - Vb**2
check("condition is w^2 + S*R*w - Vb^2 = 0 with S = u/R - u'",
      eq_A_w + target)

disc = sp.sqrt(S**2 * R_s**2 / 4 + Vb**2)
w_pro = -S * R_s / 2 + disc
w_ret = -S * R_s / 2 - disc
v_pro_A = U + w_pro
v_ret_A = U + w_ret
print(f"[A3] v_pro = u - S R/2 + sqrt(S^2 R^2/4 + Vb^2)")
print(f"     v_ret = u - S R/2 - sqrt(S^2 R^2/4 + Vb^2)")

# Solid-body flow u = Omega R has zero shear: additivity is EXACT.
check("solid-body flow (S=0): v_pro = Vb + Omega R exactly",
      (v_pro_A.subs({S: 0, U: Om * R_s})) - (Vb + Om * R_s))
check("no medium (u=0, S=0): v_pro = Vb", v_pro_A.subs({S: 0, U: 0}) - Vb)

# Chirality: Delta V = v_pro - |v_ret| = v_pro + v_ret (v_ret < 0 when the
# retrograde orbit is genuinely retrograde in the coordinate frame).
dV_A = sp.simplify(v_pro_A + v_ret_A)
print(f"[A4] Delta V = v_pro + v_ret = {dV_A}")
# 2u - S R = 2u - u + u' R = u + u' R = d(R u)/dR = R * zeta,
# where zeta = (1/R) d(R u)/dR is the vorticity of the flow.
zeta = (U + Up * R_s) / R_s
check("Delta V = R * zeta  (zeta = flow vorticity)",
      dV_A.subs(S, U / R_s - Up) - R_s * zeta)

# First order in the shear. The additive relation v = Vb + u is exact only
# for S = 0; in general v_pro = Vb + u - S R/2 + O(S^2). Since
# Delta V = 2u - S R, the first-order-consistent prediction from an observed
# rotation curve is Delta V = 2 (v_pro - Vb) for ANY flow profile, and the
# shortcut "u = V_obs - Vb, Delta V = d(R u)/dR" is wrong at first order.
ser1 = sp.series(v_pro_A, S, 0, 2).removeO()
ser2 = sp.series(v_pro_A, S, 0, 3).removeO()
print(f"\n[A6] v_pro = {sp.simplify(ser1)} + O(S^2)")
check("v_pro = Vb + u - S R/2 + O(S^2)", ser1 - (Vb + U - S * R_s / 2))
check("Delta V = 2 (v_pro - Vb) + O(S^2) for any flow profile",
      dV_A - 2 * (ser1 - Vb))
check("second-order term in v_pro is S^2 R^2 / (8 Vb)",
      sp.simplify(ser2 - ser1) - S**2 * R_s**2 / (8 * Vb))

# The RT profile as a flow: u = omega R / (1 + R/R_t).
u_RT = om * R_s / (1 + R_s / Rt)
zeta_RT = sp.simplify(sp.diff(R_s * u_RT, R_s) / R_s)
dV_RT = sp.simplify(R_s * zeta_RT)
print(f"\n[A5] RT flow u = omega R/(1+R/R_t):")
print(f"     vorticity zeta(R) = {sp.factor(zeta_RT)}")
print(f"     Delta V(R)        = {sp.factor(dV_RT)}")
print(f"     Delta V(R -> 0)   = {sp.limit(dV_RT / R_s, R_s, 0)} * R")
print(f"     Delta V(R -> oo)  = {sp.limit(dV_RT, R_s, sp.oo)}  (= V_sat)")
check("RT flow: Delta V -> 2 omega R at small R",
      sp.limit(dV_RT / R_s, R_s, 0) - 2 * om)
check("RT flow: Delta V -> omega R_t = V_sat at large R",
      sp.limit(dV_RT, R_s, sp.oo) - om * Rt)


print("\n" + "=" * 72)
print("Formulation B: vortical force field, L = 1/2 v^2 - Phi + v.A")
print("=" * 72)

L_B = sp.Rational(1, 2) * (Rd**2 + (Rf * phid)**2) - Phi(Rf) + Rf * phid * A(Rf)
eq_B = radial_circular_equation(L_B)
Bz = sp.Symbol("B")  # (curl A)_z = (1/R) d(R A_phi)/dR
eq_B = sp.expand(eq_B.subs(sp.Derivative(A(R_s), R_s), Bz - A(R_s) / R_s))
print(f"\n[B1] Circular-orbit condition (=0):\n     {sp.simplify(eq_B)}")
check("condition is v^2/R = Vb^2/R - B v  (B parallel to L slows orbit)",
      sp.simplify(eq_B * R_s) - sp.simplify(-(v**2 - Vb**2 + Bz * v * R_s)))

# For the force to *raise* prograde speeds, B must be anti-parallel to the
# disk's angular momentum. Write B = -2*Omega (Omega > 0 boosts prograde).
v_sols = sp.solve(sp.Eq(v**2, Vb**2 + 2 * Om * v * R_s), v)
v_pro_B = [s for s in v_sols if s.subs({Om: 1, R_s: 1, Vb: 1}) > 0][0]
print(f"[B2] With B = -2 Omega: v_pro = {v_pro_B}")
check("B: v_pro = Omega R + sqrt(Omega^2 R^2 + Vb^2)",
      v_pro_B - (Om * R_s + sp.sqrt(Om**2 * R_s**2 + Vb**2)))
series = sp.series(v_pro_B, Om, 0, 2).removeO()
check("B: weak field gives FC25 linear, v = Vb + Omega R + O(Omega^2)",
      series - (Vb + Om * R_s))
# Retrograde orbit: velocity reverses, so the force sign flips.
v_ret_B = -Om * R_s + sp.sqrt(Om**2 * R_s**2 + Vb**2)  # |v_ret|
dV_B = sp.simplify(v_pro_B - v_ret_B)
print(f"[B3] Delta V = v_pro - |v_ret| = {dV_B}")
check("B: Delta V = 2 Omega R = R * |B|", dV_B - 2 * Om * R_s)

print("\n" + "=" * 72)
print("Formulation A with pressure support: chirality estimator for real disks")
print("=" * 72)
# Stellar disks are not cold: their mean speed lags the circular speed by
# asymmetric drift. In the medium frame (solid-body flow u = Omega R, so
# S = 0) the Jeans equation reads  w^2 + k sigma^2 = Vb^2  for each disk,
# with k = -(dln nu/dlnR + dln sigma_R^2/dlnR + 1 - sigma_phi^2/sigma_R^2)
# * (sigma_R^2 / sigma_los^2), a geometric factor of order 1-3.
V1, V2, s1, s2, kk, uu = sp.symbols("V_1 V_2 sigma_1 sigma_2 k u",
                                    positive=True)
# Disk 1 observed at +V1, disk 2 at -V2 (counter-rotating). Speeds relative
# to the medium are V1 - u and V2 + u.
jeans_balance = sp.Eq((V1 - uu)**2 + kk * s1**2, (V2 + uu)**2 + kk * s2**2)
u_sol = sp.solve(jeans_balance, uu)
check("Jeans balance is linear in u (unique solution)", len(u_sol) - 1)
u_hat = u_sol[0]
dV_hat = sp.simplify(2 * u_hat)
print(f"\n[C1] R*zeta = 2u = {dV_hat}")
check("estimator: R*zeta = V1 - V2 - k (s2^2 - s1^2)/(V1 + V2)",
      dV_hat - (V1 - V2 - kk * (s2**2 - s1**2) / (V1 + V2)))
check("cold disks (sigma=0) reduce to Delta V = V1 - V2",
      dV_hat.subs({s1: 0, s2: 0}) - (V1 - V2))
# Mean speed relative to the medium recovers the baryonic circular speed:
mean_rel = sp.simplify(((V1 - u_hat) + (V2 + u_hat)) / 2)
check("mean of medium-frame speeds = (V1 + V2)/2 (u cancels)",
      mean_rel - (V1 + V2) / 2)
print("     => in a counter-rotating galaxy the MEAN of the two (drift-")
print("        corrected) speeds measures V_bary, and the whole 'missing")
print("        mass' signal appears as the asymmetry R*zeta.")

print("\n" + "=" * 72)
print("Summary")
print("=" * 72)
print("""
 * Both formulations predict a prograde/retrograde speed asymmetry
       Delta V(R) = R * zeta(R)
   where zeta is the local vorticity of space (A) or |curl A| (B).
   Dark matter and MOND predict Delta V = 0 (after asymmetric-drift
   correction), because their forces depend only on position.

 * Formulation A makes additive coupling EXACT for solid-body flow, and
   identifies the RT correction with the flow speed of space, u(R).
   For the RT profile, Delta V rises from 2*omega*R to V_sat = omega*R_t.
   For a general flow, v_pro = Vb + u - S R/2 + O(S^2), so the prediction
   from an observed rotation curve is Delta V = 2 (V_obs - Vb) + O(S^2)
   whatever the profile (the same as Formulation B and as rigid dragging).

 * Neither formulation adds a radial force on an isotropic (pressure-
   supported) orbit population, and neither adds lensing mass.
   (Gate G2.)
""")
print("All assertions passed.")
