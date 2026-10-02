# Gate G1 with MaNGA data (REGISTERED 2026-10-02 as PREREG amendment c; this file is the review copy)

This draft is for review. Once agreed, it will be appended to PREREG.md as a dated amendment, *before any MaNGA velocity amplitudes are looked at*.

## Question

Do galaxies whose gas counter-rotates against their stars show a speed asymmetry R·ζ that matched co-rotating galaxies do not?

## Data (public, SDSS DR17)

- **MaNGA DAP MAPS files.** Type HYB10-MILESHC-MASTARSSP, about 10–30 MB per galaxy. They contain stellar V and σ, and Hα gas V and σ, per spaxel, with errors and quality masks.
- **Catalogues:**
  - `drpall` and `dapall` (redshift, R_e, NSA ellipticity, stellar mass);
  - MaNGA PyMorph and Deep-Learning morphology value-added catalogues, for T-type and disk scale length.

## Sample selection (uses position angles only; amplitudes stay unseen)

1. **Quality.**
   - At least 50 valid spaxels in both the stellar and Hα velocity maps (DAP masks, S/N ≥ 5).
   - Photometric inclination between 30° and 80° (NSA b/a, with intrinsic thickness q₀ = 0.2).
   - Not flagged as a merger or interacting galaxy in the morphology catalogue.
2. **Kinematic position angles.** PA_star and PA_gas come from `fit_kinematic_pa` (Krajnović+2006). ΔPA = |PA_gas − PA_star|, folded to 0–180°.
3. **Groups.**
   - Counter-rotators (CR): ΔPA > 150°.
   - Co-rotators (CO): ΔPA < 30°.
   - Galaxies in between are dropped (misaligned, not counter-rotating).
4. **Matching.** For each CR galaxy, pick the 3 nearest CO galaxies in (log M*, T-type, log σ_*,Re, inclination). Matching is done before any rotation amplitude is computed.

## Measurement (identical pipeline for both groups)

1. **Rotation curves.** Fit tilted-ring (or harmonic) curves along the common kinematic PA to the stellar and gas velocity fields. Deproject with the photometric inclination.
2. **Radii.** Evaluate at R/R_e = 0.5, 1.0 and 1.5, using only radii where both components have valid rings.
3. **Gas dispersion.** Correct σ_gas for instrumental broadening (DAP provides it).
4. **Asymmetric drift** for both components:

   k = R/h_R + R/h_σ − ½

   Here h_R is the photometric disk scale length and h_σ is the dispersion scale length fitted to σ(R). The systematic range is k ± 1.
5. **Per-galaxy statistic at each radius:**

   D = (V_gas² + k_g σ_gas² − V_*² − k_* σ_*²) / (V_gas + V_*)

   - D is the drift-corrected circular speed of the gas minus that of the stars.
   - For co-rotators, H_vort predicts D ≈ 0, since both components see the same flow.
   - For counter-rotators, H_vort predicts |D| ≈ R·ζ.
   - Gas disequilibrium and M/L systematics affect both groups and cancel in the comparison.

## Tests

**Primary test (sign-agnostic).**
- Excess variance: E = Var(D_CR) − Var(D_CO).
- H_vort predicts E ≈ ⟨(R·ζ)²⟩. The expected value comes from the SPARC RT-flow vorticity at the same R and V, as in the literature version of G1.
- Uncertainty: bootstrap over galaxies (10,000 resamples).

**Secondary test (signed).** Assume the flow co-rotates with the dominant stellar disk. Then D_CR < 0, i.e. the gas runs slow by R·ζ. Statistic: median(D_CR) − median(D_CO).

**Verdicts at R = 1 R_e** (pre-specified primary radius):
- *Falsified at the predicted amplitude:* the upper 2σ bound on E is below 25% of the predicted E.
- *Supported:* E > 0 at more than 3σ and consistent with the prediction within 2σ.
- *Otherwise:* inconclusive. Report the measured E and the minimum sample size needed.

## Power check (to run on synthetic data before any real data)

Inject R·ζ drawn from the SPARC expectation into mock CR galaxies, using realistic MaNGA errors and k systematics. Report the number of CR galaxies needed for a 3σ detection.

Published samples (Jin+2016; Chen+2016) suggest about 50–100 gas–star counter-rotators in DR17.

## Decisions for the reviewer

1. **Primary statistic.** Is the sign-agnostic variance test the right primary? The signed test is more powerful but assumes the flow follows the stars.
2. **ΔPA cut.** Is 150° acceptable? Bevacqua+2022 and ATLAS³D use similar cuts.
3. **Bias from slow early-type galaxies.** Low-V galaxies give small R·ζ predictions. Should the sample be restricted to galaxies with V_* > 100 km/s, where the predicted signal is largest?
