#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Order-of-magnitude gates for the three hypotheses.

Each gate compares what a mechanism can deliver with what the rotation
curves require: a vorticity / angular rate of order omega ~ 10-50 km/s/kpc
at R ~ 1-10 kpc (SPARC median omega for taper-preferred galaxies is of this
order; see notebooks/00_foundation.ipynb for the measured distribution).

  G-seed : GR frame-dragging from a central SMBH at kpc radii.
  G-spin : Primordial cosmic vorticity (CMB bound), amplified by collapse
           under Kelvin's circulation theorem.
  G-vort : (a) the coupling a current-sourced vorticity field needs;
           (b) the solar-system frame-dragging that coupling would imply;
           (c) the a0/V magnitude coincidence;
           (d) Gate G3 -- rotation of local inertial frames at the Sun
               implied by the Milky Way's own flow, vs Lunar Laser Ranging.

All inputs are listed with sources at the top. Exit code 0 = script ran and
every gate printed a verdict (the verdicts themselves are physics results,
not assertions).
"""

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Constants (SI) and unit conversions
# ---------------------------------------------------------------------------
G = 6.674e-11            # m^3 kg^-1 s^-2
c = 2.998e8              # m/s
M_SUN = 1.989e30         # kg
PC = 3.0857e16           # m
KPC = 1e3 * PC
MPC = 1e6 * PC
YR = 3.156e7             # s
KMS_PER_KPC = 1e3 / KPC  # 1 km/s/kpc in s^-1 (= 3.24e-17)
MAS = np.pi / (180 * 3600 * 1e3)  # 1 mas in rad
A0 = 1.2e-10             # m/s^2, MOND acceleration scale
H0 = 70e3 / MPC          # s^-1

def to_kms_kpc(rate_si):
    return rate_si / KMS_PER_KPC

def to_mas_yr(rate_si):
    return rate_si * YR / MAS

# Requirement from rotation curves (order of magnitude).
OMEGA_REQ = 20.0 * KMS_PER_KPC   # s^-1

def verdict(name, delivered, required):
    ratio = delivered / required
    status = "PASSES" if ratio >= 0.1 else "FAILS"
    print(f"\n    >>> {name}: delivers {ratio:.1e} of the requirement"
          f" -> {status} (gap {1/ratio:.1e}x)" if ratio < 1 else
          f"\n    >>> {name}: delivers {ratio:.1e} of the requirement -> PASSES")
    return ratio


print("=" * 72)
print(f"Requirement: omega ~ 20 km/s/kpc = {OMEGA_REQ:.2e} s^-1")
print("=" * 72)

# ---------------------------------------------------------------------------
# G-seed: Lense-Thirring from a maximally spinning SMBH
# ---------------------------------------------------------------------------
print("\n[G-seed] Frame-dragging from a central SMBH (GR)")
for M_bh in (1e6, 1e8, 1e10):
    M = M_bh * M_SUN
    J = 1.0 * G * M**2 / c          # spin parameter a* = 1
    r = 1 * KPC
    om_lt = 2 * G * J / (c**2 * r**3)
    print(f"    M_BH = 1e{int(np.log10(M_bh))} Msun, r = 1 kpc:"
          f" Omega_LT = {om_lt:.1e} s^-1 = {to_kms_kpc(om_lt):.1e} km/s/kpc")
M = 1e10 * M_SUN
verdict("SMBH frame-dragging (1e10 Msun, a*=1, 1 kpc)",
        2 * G * (G * M**2 / c) / (c**2 * KPC**3), OMEGA_REQ)
print("    H_seed cannot supply the amplitude; it can only shape structure.")

# ---------------------------------------------------------------------------
# G-spin: primordial cosmic vorticity
# ---------------------------------------------------------------------------
print("\n[G-spin] Primordial cosmic vorticity amplified by collapse")
OMEGA_OVER_H_BOUND = 4.7e-11     # Saadeh et al. 2016 (PRL 117, 131302), 95%
zeta0 = OMEGA_OVER_H_BOUND * H0
print(f"    CMB bound: (omega/H)_0 < {OMEGA_OVER_H_BOUND:.1e}"
      f" -> zeta_0 < {zeta0:.1e} s^-1")
# Kelvin: circulation conserved, zeta * area = const. Physical vorticity of
# the background scales as a^-2 and so does 1/area of a comoving patch, so
# the redshift of collapse cancels; only the comoving Lagrangian radius
# matters.
rho_m = 0.3 * 3 * H0**2 / (8 * np.pi * G)      # mean matter density today
for M_halo in (1e10, 1e12):
    R_L = (3 * M_halo * M_SUN / (4 * np.pi * rho_m)) ** (1 / 3)
    R_disk = 3 * KPC
    amp = (R_L / R_disk) ** 2
    print(f"    M_halo = 1e{int(np.log10(M_halo))}: Lagrangian radius"
          f" {R_L / MPC:.2f} Mpc, collapse to 3 kpc amplifies by {amp:.1e}"
          f" -> zeta_disk < {zeta0 * amp:.1e} s^-1")
R_L = (3 * 1e12 * M_SUN / (4 * np.pi * rho_m)) ** (1 / 3)
verdict("Primordial vorticity after Kelvin amplification",
        zeta0 * (R_L / (3 * KPC)) ** 2, OMEGA_REQ)
print("    H_spin cannot supply the amplitude; at most a handedness bias.")

# ---------------------------------------------------------------------------
# G-vort (a): coupling needed for a current-sourced field
# ---------------------------------------------------------------------------
print("\n[G-vort a] Coupling for curl(Omega) = kappa * rho * v")
Sigma = 100 * M_SUN / PC**2     # typical disk surface density
v_rot = 150e3
# Thin current sheet: field jump = kappa * Sigma * v, so Omega ~ kappa Sigma v / 2
kappa_req = 2 * OMEGA_REQ / (Sigma * v_rot)
kappa_gr = G / c**2
c_v = np.sqrt(G / kappa_req)
print(f"    Sigma = 100 Msun/pc^2, v = 150 km/s")
print(f"    kappa required = {kappa_req:.1e} m/kg;  G/c^2 = {kappa_gr:.1e} m/kg")
print(f"    enhancement over GR = {kappa_req / kappa_gr:.1e}"
      f"  (MGEM measured eta ~ 2e8)")
print(f"    equivalent characteristic speed c_v = sqrt(G/kappa) ="
      f" {c_v / 1e3:.0f} km/s")

# ---------------------------------------------------------------------------
# G-vort (b): what that coupling does to Earth's frame-dragging
# ---------------------------------------------------------------------------
print("\n[G-vort b] Earth frame-dragging with the same coupling")
GPB_GR = 39.2       # mas/yr, GR prediction (Everitt et al. 2011)
GPB_ERR = 7.2       # mas/yr, measurement uncertainty
enh = kappa_req / kappa_gr
print(f"    GR: {GPB_GR} mas/yr; measured to +/-{GPB_ERR} mas/yr (Gravity Probe B)")
print(f"    Unscreened enhanced coupling predicts {GPB_GR * enh:.1e} mas/yr")
print(f"    >>> Excluded by {GPB_GR * enh / GPB_ERR:.0e} sigma unless the field"
      f" is screened at high acceleration (g_Earth ~ {G*5.97e24/(7.0e6)**2:.0f}"
      f" m/s^2 >> a0).")

# ---------------------------------------------------------------------------
# G-vort (c): the a0/V coincidence
# ---------------------------------------------------------------------------
print("\n[G-vort c] Magnitude of Omega ~ a0 / V")
for V in (50, 100, 150, 250):
    print(f"    V = {V:3d} km/s: a0/V = {to_kms_kpc(A0 / (V * 1e3)):5.1f} km/s/kpc")
print(f"    (a0 / (c H0) = {A0 / (c * H0):.3f}; 1/(2 pi) = {1/(2*np.pi):.3f})")
print("    The only constant that reproduces omega is a0 -> any viable vortex")
print("    theory carries an acceleration scale, and screening at g >> a0 is the")
print("    natural escape from [G-vort b] and [G3].")

# ---------------------------------------------------------------------------
# G3: Milky Way flow at the Sun vs Lunar Laser Ranging
# ---------------------------------------------------------------------------
print("\n[G3] Local inertial-frame rotation at the Sun (unscreened H_vort)")
R0 = 8.122                       # kpc, GRAVITY 2018
VC = 229.0                       # km/s, Eilers et al. 2019
LLR_TIE_RATE = 0.02              # mas/yr, dynamical-ICRF tie rate (Hofmann+)
LLR_GP_ERR = 0.12                # mas/yr, geodetic precession uncertainty
print(f"    V_c(R0 = {R0} kpc) = {VC} km/s")
print(f"    Bounds: dynamical-vs-quasar frame drift {LLR_TIE_RATE} mas/yr;"
      f" geodetic precession +/-{LLR_GP_ERR} mas/yr")
print(f"    {'V_bary':>7} {'u=Vc-Vb':>8} {'zeta/2 [km/s/kpc]':>18}"
      f" {'precession [mas/yr]':>20} {'x LLR GP bound':>15}")
precessions = []
for Vb in (150, 175, 200):
    u = VC - Vb
    # Fluid elements rotate at zeta/2. zeta = (1/R) d(R u)/dR lies between
    # u/R (flat u) and 2u/R (solid body).
    lo, hi = u / R0 / 2, u / R0
    p_lo = to_mas_yr(lo * KMS_PER_KPC)
    p_hi = to_mas_yr(hi * KMS_PER_KPC)
    precessions += [p_lo, p_hi]
    print(f"    {Vb:7d} {u:8.0f} {lo:8.1f} - {hi:<8.1f}"
          f" {p_lo:9.2f} - {p_hi:<9.2f} {p_lo / LLR_GP_ERR:6.0f} - {p_hi / LLR_GP_ERR:<6.0f}")
print(f"    >>> Unscreened H_vort predicts {min(precessions):.2f}-"
      f"{max(precessions):.2f} mas/yr; LLR allows <~{LLR_GP_ERR}.")
print("        Excluded unless screened inside the solar system, where")
print(f"        g_sun(1 AU) = {G * M_SUN / (1.496e11)**2:.1e} m/s^2 = "
      f"{G * M_SUN / (1.496e11)**2 / A0:.0e} a0.")

# ---------------------------------------------------------------------------
# Dark energy inside galaxies
# ---------------------------------------------------------------------------
print("\n[Lambda] Cosmological-constant acceleration inside a galaxy")
OMEGA_L = 0.7
# Newtonian limit of GR with Lambda: a_Lambda = +(Lambda c^2 / 3) r
#                                             = Omega_L H0^2 r  (outward, radial)
for r_kpc in (1, 10, 30, 100):
    aL = OMEGA_L * H0**2 * r_kpc * KPC
    print(f"    r = {r_kpc:4d} kpc: a_Lambda = {aL:.1e} m/s^2 = {aL / A0:.1e} a0")
# Radius where Lambda repulsion equals the attraction of mass M
for M in (1e10, 1e12):
    r_zg = (G * M * M_SUN / (OMEGA_L * H0**2)) ** (1 / 3)
    print(f"    zero-gravity radius for M = 1e{int(np.log10(M))} Msun:"
          f" {r_zg / MPC:.2f} Mpc")
print(f"    Coincidence: c * sqrt(Lambda/3) / (2 pi) = c H0 sqrt(Omega_L)/(2 pi)"
      f" = {c * H0 * np.sqrt(OMEGA_L) / (2 * np.pi):.1e} m/s^2 vs a0 = {A0:.1e}")
print("    Lambda is ~1e-5 of a0 at 10 kpc, radial and velocity-independent:")
print("    no prograde/retrograde asymmetry and no role in disk dynamics. It")
print("    matters only at the ~Mpc zero-gravity radius (group outskirts).")

print("\nDone.")
