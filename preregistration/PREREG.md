# Pre-registration: vortical-vacuum and seed-first tests

Registered: 2026-10-02, before running any test below on data.
Author: J. Schneider (with Claude Code). File is append-only: later changes go in a dated "Amendments" section at the bottom, never as edits above it.

**Disclosure.** Before registration, published values for NGC 4550 were viewed while checking that Gate G1 is feasible (Johnston+2013; Coccato+2013). Results for NGC 4550 are therefore reported but flagged as *not blind*.

## Common definitions

- **Clean sample.** SPARC galaxies with:
  - Q < 3;
  - Inc ≥ 30°;
  - RT fit flag OK;
  - R_t identifiable: σ(R_t) < R_t and R_t < R_max, where R_max is the last measured radius.
- **Velocities** are always fit to V_obs with SPARC errors. Υ_disk = 0.5 and Υ_bulge = 0.7, both fixed.
- **Units.** Π₂ = G M_bar / (a₀ R_d²), with M_bar = 1.33 M_HI + Υ_d L_disk + Υ_b L_bul and a₀ = 1.2×10⁻¹⁰ m s⁻².
- **Significance.** All correlation p-values come from 10,000 permutations (shuffle nulls). The headline threshold is p < 0.01.

## Step 0: Foundation (data cleaning, not a hypothesis test)

Recompute on the clean sample, and also on the "taper-preferred" sample (ΔBIC > 2):
- the R_t–R_d log-log slope and its correct OLS standard error;
- the median R_t/R_d;
- a shuffle null for the R_t–R_d correlation.

## Gate G1: Counter-rotation chirality (H_vort)

**Prediction under H_vort.** Two counter-rotating disks at the same radius R have drift-corrected speeds that differ by R·ζ:

R·ζ̂ = V₁ − V₂ − k(σ₂² − σ₁²)/(V₁ + V₂)

The expected magnitude for a system is the distribution of R·ζ_RT(R) over clean SPARC galaxies with |V_flat − V_sys| < 30 km/s, evaluated at the same R. ζ_RT is the vorticity of the fitted RT flow.

**Kill criterion, per system.** H_vort, as the sole source of the mass discrepancy, is disfavoured for a system if the 2σ upper bound of |R·ζ̂| is below the 16th percentile of the expected distribution. The bound is taken as the maximum over k ∈ {1, 2, 3}.

**Kill criterion, overall.** H_vort is falsified as the sole source if this holds for at least 3 systems, or for every system tested when fewer than 3 are available. NGC 4550 is reported but does not count towards the 3, because it is not blind.

**Known limitation.** Most counter-rotators are S0/Sa galaxies. SPARC contains few of these, so the expected distribution may overstate ζ in baryon-dominated inner regions.

## Gate G2: Pressure-supported systems (literature)

A v×Ω force, or a flowing-space frame, adds no radial support to an isotropic orbit population. H_vort as the sole source therefore predicts no mass discrepancy in dSphs or ellipticals.

This prediction is assessed against published results only (Wolf+2010; Lelli+2017; Brouwer+2021 weak-lensing RAR). **Verdict rule.** If published discrepancies in pressure-supported systems are significant at more than 5σ, H_vort is downgraded to "at most a partial contributor."

## Gate G3: Local inertial frame (calculation)

**Unscreened prediction.** Inertial frames precess at about ζ_MW/2 at the Sun.

**Bound.** Lunar Laser Ranging: geodetic precession agrees with GR to ±0.12 mas/yr, and the dynamical–ICRF tie rate is 0.02 mas/yr.

**Verdict rule.** If the predicted precession is more than 3× the bound, unscreened H_vort is excluded and any surviving version must include screening at g ≫ a₀.

## Gate G4: Flowing-space fit (Formulation A)

**Model.** V = u − SR/2 + √(S²R²/4 + V_b²), with u = ωR/(1+R/R_t). It has 2 free parameters, the same as RT-additive.

**Comparison.** Per galaxy on the clean sample: ΔBIC against RT-additive (refit with the same code), and against the NFW and MOND-free BIC from the 2026b tournament table.

**Report.** Win fractions and median ΔBIC. Formulation A is "competitive" if median ΔBIC(A − NFW) ≤ 2. There is no kill criterion: a better fit cannot validate H_vort, because the shapes are flexible. Only G1 can.

## S2: Scaling law R_t/R_d = F(Π₂, Π₃)

- **S2a (pure-a₀ rejection).** A theory with no disk length predicts ω·V_sat/a₀ = const and R_t ∝ M_bar^0.5. Rejected if the clean-sample slope of log R_t on log M_bar differs from 0.5 at more than 3σ.
- **S2b (F depends on Π₂).** Spearman ρ(R_t/R_d, Π₂) on the clean sample. Supported if p < 0.01 under permutation.
- **S2c (seed term Π₃; pending a direct M_BH cross-match).** Partial Spearman ρ(R_t/R_d, M_BH/M_bar | Π₂). H_seed is supported if p < 0.01 and the sign is positive (following the 2026d B/T result).
  - Only dynamical M_BH measurements are used: maser, stellar, gas dynamics or reverberation mapping.
  - If N < 15, report it as exploratory.

## P1: Spin-handedness dipole (H_spin; pending catalogue choice)

- **Model.** The probability that a galaxy's spin is "Z-wise" is ½(1 + A·n̂·d̂). Fit (A, d̂) by maximum likelihood. There is no grid search over directions.
- **Significance.** 10,000 sign-shuffled mocks with the same sky positions. Threshold p < 0.003.
- **Catalogue.** Fixed in an amendment before the data are downloaded.

---

## Amendments

### 2026-10-02 (after the preliminary, non-blind NGC 4550 run)

The NGC 4550 result is not decisive. With two stellar disks, the uncertainty in the asymmetric-drift factor k ∈ [1, 3] dominates. The blind G1 sample therefore adds a second estimator for **cold gas** (k ≈ 0):

- **G1-gas.** Stars and gas counter-rotate at the same R. The stellar speed is drift-corrected with k ∈ {1, 2, 3}; the gas speed is used directly. The same estimator applies, with σ_gas set to 0.
- **G1-reversal.** Galaxies whose gas reverses rotation sense with radius (e.g. NGC 4826) offer a further test. H_vort predicts that |V_gas| jumps by R·ζ across the reversal radius, because V_bary and the flow are continuous there. Estimator: the difference between |V| fitted just inside and just outside the reversal, compared with the same expected distribution. The kill criterion is unchanged.
- **Candidates**, to be fixed before their kinematics are extracted:
  - NGC 4826 (Braun+1994; Rubin 1994)
  - NGC 4138 and NGC 3626 (Jore+1996; Ciri+1995)
  - NGC 3593 (Coccato+2013). The stellar values were seen, but not used in any test.
  - ATLAS³D early-type galaxies with kinematically misaligned (≈180°) gas: Davis+2011; Davis+2013 (CO versus JAM circular speed).

### 2026-10-02 (b): G1 extraction protocol, fixed before any candidate kinematics are viewed

1. **Source rule.** For each candidate, use the most recent peer-reviewed study that gives both components' rotation speeds in text or tables. Values read off figures are deferred: the system is listed as "pending digitisation" and not tested. If several studies qualify, use the one with the largest radial coverage. Report the others as cross-checks only.
2. **Radius rule.**
   - Same-radius systems: use the outermost radius at which both components are measured, where the expected R·ζ is largest.
   - Reversal systems: interpolate each side's |V(R)| to the reversal radius R_rev, using only points within 0.3 R_rev of it.
3. **Deprojection.** Use the source's own inclination. If it gives none, use the HyperLeda value.
4. **Gas.** Treat ionised or molecular gas with σ ≤ 20 km/s as cold, with σ set to 0 in the estimator. Gas with σ > 20 km/s gets the same k range as the stars.
5. **Uncertainties.** Use the quoted errors. Where none are quoted, assume ±10 km/s for velocities and ±25% for dispersions.
6. **Reversal estimator.** R·ζ̂ = |V_out(R_rev)| − |V_in(R_rev)|. The sign is not predicted, because the flow direction is unknown. The kill criterion applies to |R·ζ̂|.
7. **Exclusions.** The only allowed reasons are named in advance: the source flags non-circular motion at the chosen radius (bar, warp, or outflow), or the components lie in different planes by more than 30°. Any exclusion is reported, with its reason.

### 2026-10-02 (c): Gate G1 with MaNGA integral-field data. REGISTERED before any MaNGA data are downloaded

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

**Reviewer decisions (J. Schneider, 2026-10-02):**
1. **Primary test:** sign-agnostic excess variance E. The signed test is secondary.
2. **Speed cut:** V_* > 100 km/s at 1 R_e. The cut is applied automatically in code, identically to the counter-rotating and co-rotating groups, before D is computed. No amplitudes are inspected by hand before the test runs.
3. **Counter-rotator definition:** ΔPA > 150°.

The power check on synthetic data is run, and its result recorded here, before any real MaNGA kinematics are loaded.

### 2026-10-02 (d): MaNGA G1 power-check result and disturbance control (before any MaNGA data)

**Power check** (`analysis/g1_manga_power.py`, synthetic data only).
- Expected R·ζ at R_eff for SPARC galaxies with V_flat > 100: median 61 km/s, RMS 69 km/s.
- Predicted E = 4730 km²/s². Shared nuisance per galaxy ≈ 30 km/s.
- If counter-rotating gas is no more disturbed than co-rotating gas, N_CR = 40 gives:
  - 3σ detection probability 1.00 when the signal is present;
  - falsification probability 0.99 when it is absent;
  - false-positive rate 0.00.
- If counter-rotating gas is 1.3× more disturbed, the 3σ false-positive rate for "E > 0" rises with N: 0.03 at N = 40, 0.21 at N = 100, 0.50 at N = 150.
- In that case the "supported" verdict is still protected by the registered requirement that E agree with the prediction within 2σ. Disturbance alone gives E ≈ 500 against 4730 predicted.

**Added control (registered now).**
- For every galaxy, compute a gas-disturbance index that does not depend on rotation amplitude: the RMS of the gas velocity-field residuals after the tilted-ring model, divided by the median gas velocity error, plus the radial twist of the gas kinematic PA.
- Report E in the full sample, and also restricted to galaxies below the combined median of the disturbance index.
- If E is significant only in the full sample, and not in the low-disturbance half at a sample size where the power check predicts detection, the result is classified **"disturbance-driven, not H_vort"**.

**Minimum sample.** N_CR ≥ 40 after all cuts. Below that, the result is reported as underpowered and inconclusive.

### 2026-10-02 (e): MaNGA data route and selection implementation (before any maps are fetched)

**Data route.** Per-galaxy DAP maps (DR17, HYB10-MILESHC-MASTARSSP) are fetched from the SDSS Marvin API with verified TLS, at about 25 KB per map. This makes it affordable to fit ΔPA for the **whole** parent sample ourselves, exactly as amendment (c) specifies. No external misalignment catalogue is used.

**Parent sample.** Built from drpall v3_1_1 and dapall 3.1.0:
- Main-sample targets only: MNGTARG1 bits 10, 11 or 12.
- DRP3QUAL and DAPQUAL with no critical bit set.
- One observation per MaNGA ID.
- Inclination between 30° and 80°, from NSA elliptical-Petrosian b/a with q₀ = 0.2.
- This gives N = 7582.

**Valid spaxels.** A spaxel counts as valid if:
- its DAP mask is 0 and its ivar is > 0;
- its velocity error is ≤ 30 km/s;
- for gas only, additionally, Hα flux S/N ≥ 5.

**Kinematic PA.**
- Use `pafit.fit_kinematic_pa` (Krajnović+2006). Coordinates are spaxel offsets from the map centre, in arcsec. The median velocity is subtracted first.
- pafit returns 0–180°. The receding side is fixed by adding 180° when Σ v·(projected distance along the PA) < 0.
- ΔPA = |PA_gas − PA_star|, folded to 0–180°. The same convention for both components cancels any axis-orientation choice.

**Groups.** Counter-rotators: ΔPA > 150°. Co-rotators: ΔPA < 30°. A galaxy is used only if both PA errors are ≤ 20°.

**Outputs at this stage.** PAs, PA errors, spaxel counts and ΔPA only. No rotation amplitudes are computed or written until the measurement stage, which applies the V_* > 100 km/s cut automatically.

*Implementation note to (e), 2026-10-02.* The PA fits use every second spaxel in x and y. Adjacent spaxels are correlated (the PSF is about 5 spaxels wide), and this makes the fits about 4× faster. The 50-valid-spaxel minimum is counted after this decimation, so it is stricter. A 20-galaxy pipeline test produced PAs only, with no amplitudes; those outputs are discarded and recomputed.

### 2026-10-02 (f): MaNGA G1 stage-2 measurement details (registered before any amplitude is measured)

1. **Merger and interaction exclusion.** Use the Galaxy Zoo MaNGA catalogue (MaNGA_gz-v2_0_1): exclude a galaxy if `t08_odd_feature_a24_merger_flag` = 1 or `t08_odd_feature_a21_disturbed_flag` = 1. T-type comes from the deep-learning morphology catalogue (manga-morphology-dl-DR17).
2. **Rings.**
   - Elliptical annuli, using the photometric inclination and each component's own kinematic PA.
   - Width 0.25 R_e, centred at 0.5, 0.75, 1.0, 1.25 and 1.5 R_e.
   - Per ring, fit v = V_sys + c₁cos θ + s₁sin θ to all valid (non-decimated) spaxels, then V_rot = c₁/sin i.
   - A ring is valid only if it has ≥ 20 spaxels and coverage in ≥ 3 of 4 azimuthal quadrants.
3. **Dispersions.**
   - Stars: √(σ² − σ_corr²), using the DAP `stellar_sigmacorr`.
   - Gas: √(σ² − σ_inst²), using the DAP `emline_instsigma`.
   - Ring value: the median.
4. **Asymmetric-drift factor: clarification of (c).** For each component,

   k = R/h_R + R/h_{σ²} − ½

   - h_R = R_e/1.678 (exponential disk).
   - h_{σ²} is the e-folding length of σ²(R), fitted to the five rings (this is the exact Jeans form for a flat curve).
   - The registered wording "h_σ" is read as h_{σ²}. Using the e-folding length of σ itself would drop a factor of 2.
   - k is clipped to the range [0, 4].
   - The systematic range k ± 1 is applied to both components together. Every verdict must hold at k − 1, k and k + 1.
5. **Speed cut.** Stellar V_rot(1 R_e) > 100 km/s, applied in code to both groups before matching.
6. **Matching.**
   - Standardised Euclidean distance on (log M*, T-type, log σ_*,1Re from dapall `STELLAR_SIGMA_1RE`, inclination).
   - Greedy, without replacement, 3 controls per counter-rotator.
   - Counter-rotators are processed in a random order (seed 20261002).
7. **Prediction E_pred.** For each counter-rotator, take the mean of (R·ζ_RT)² over clean SPARC galaxies with |V_flat − V_c,*| < 30 km/s, evaluated at R = R_e in kpc. Here V_c,* is the galaxy's drift-corrected stellar speed at 1 R_e. E_pred is the mean of these values over all counter-rotators.
8. **Disturbance index.** The mean of the percentile ranks of two quantities:
   - (a) RMS gas residual about the ring models (0.5–1.5 R_e), divided by the median gas velocity error;
   - (b) the circular difference between the gas PA fitted to spaxels inside R_e and the gas PA fitted to spaxels outside it.
9. **Pipeline validation.** Before running on real data, the pipeline must recover an injected R·ζ from synthetic velocity and dispersion maps to within 10%.

*Validation note to (f).9, 2026-10-02: PASSED on synthetic data* (`tests/test_manga_measure.py`). The mock is a counter-rotating gas and star system in a solid-body flow.
- At 1 R_e the pipeline recovers D = −34.7 against a truth of −36.0 km/s, and −71.3 against −72.0.
- With no flow, |D| ≤ 2 km/s at every ring.
- A shared k shift of ±1 moves D by about ±8–12 km/s.

*Data-route note to (e) and (f).3, 2026-10-02.*
- **Route change.** The Marvin API rate-limited the run (HTTP 429) after about 140 galaxies, so maps are now extracted from the official SAS MAPS files (data.sdss.org, same DAP version and type). The extracted stellar velocity, Hα velocity and Hα flux maps (value, ivar and mask) were verified identical to the Marvin data for a test galaxy.
- **σ-correction channel.** STELLAR_SIGMACORR uses channel C1 ('resolution difference'), following the DR17 DAP guidance (sdss4.org/dr17/manga/manga-data/working-with-manga-data).

### 2026-10-02 (g): MaNGA G1 final design (registered before any D is computed; supersedes the primary test and speed cut of (c))

**Why the change.** A realistic power check (`analysis/g1_manga_power_real.py` and `analysis/g1_manga_power_signed.py`) used the real sample's nuisance properties: stellar V and σ, gas σ, k, ring errors and disturbance. It did not use gas rotation speeds or D.
- Only 20 counter-rotators pass V_* > 100 km/s.
- The registered sign-agnostic variance test has about 0% probability of reaching any verdict at cuts from 40 to 100 km/s.
- The signed test does far better (table below).

Counts were seen before this amendment. D and E were not.

**Physical justification for the sign.** Counter-rotating gas is a minority component accreted later. A flow sourced by mass currents, or inherited from the galaxy's formation (H_seed), follows the stellar disk, which carries nearly all the angular momentum. The signed test therefore assumes the flow co-rotates with the stars. Counter-rotating gas is then predicted to run slow by R·ζ.

**Primary test.**
- Statistic: S = median(D_CO) − median(D_CR) at 1 R_e, with 10,000 bootstrap resamples.
- Prediction: S_pred = mean over counter-rotators of the mean R·ζ_RT, taken from clean SPARC galaxies with |V_flat − V_c,*| < 30 km/s and evaluated at R_e in kpc.
- Sample cut: stellar V_rot(1 R_e) > **40 km/s**, applied to both groups before matching. This gives N_CR = 81 and N_CO = 243.

**Verdicts.** Each must hold at k − 1, k and k + 1.
- **FALSIFIED** (a flow following the stars, at ≥ ½ the RT-implied amplitude): S + 2σ < **0.5 S_pred**.
- **SUPPORTED:** S > 3σ and |S − S_pred| < 2σ.
- **DISTURBANCE-DRIVEN:** S > 3σ in the full sample but not in the low-disturbance half, provided that half has ≥ 20 counter-rotators.
- **INCONCLUSIVE:** anything else. The minimum N_CR remains 40.

**Simulated power for this design.**

| Truth | Outcome |
|---|---|
| Signal present | SUPPORTED with probability 0.89 |
| Signal absent | FALSIFIED with probability 0.94 |
| Either | No wrong verdicts in 300 simulations of each case |

**Secondary tests (reported, not used for the verdict).**
- The sign-agnostic E (underpowered).
- All results at 0.5 and 1.5 R_e.
- The low-disturbance subsample.

**Scope of a FALSIFIED verdict.** It excludes a stellar-aligned flow at ≥ 50% of the amplitude implied by the RT fits. Weaker flows, or flows aligned with the gas, are not excluded.

---

## Results

### 2026-10-02: MaNGA G1 (design g), unblinded once

**Sample.** N_CR = 81 and N_CO = 243, as registered. Log: `results/g1_manga_test.log`. Output: `results/g1_manga.json`.

**Primary test at 1 R_e.** S_pred = 47.3 km/s.

| | k − 1 | k | k + 1 |
|---|---|---|---|
| S (km/s) | +13.9 ± 9.3 | +26.0 ± 9.5 | +33.4 ± 10.3 |
| Low-disturbance half (N_CR = 34) | +7.9 ± 12.1 | +13.9 ± 12.6 | +23.0 ± 10.8 |
| Median D, co-rotators (km/s) | +11.7 | −2.9 | −19.5 |

**Registered verdict: INCONCLUSIVE.**
- Not SUPPORTED: S does not exceed 3σ at all three k (1.5σ, 2.7σ, 3.2σ).
- Not FALSIFIED: S + 2σ exceeds 0.5 S_pred = 23.7 at every k.
- Not DISTURBANCE-DRIVEN: the full-sample S is not above 3σ at all k.

**Secondary results.**
- At 0.5 R_e: S = +0.3 to +35.3, with large errors.
- At 1.5 R_e: S = +29.6 ± 10.9, +28.8 ± 11.4 and +31.3 ± 15.0. In the low-disturbance half: +29.8 to +41.9, with errors of 15–17.
- The sign-agnostic E is dominated by heavy-tailed outliers (E ~ 10⁴–10⁵ km²/s²) and is uninformative.

**Notes (not part of the registered rules).**
- The observed sign matches H_vort with a stellar-aligned flow. It also matches a conventional explanation: accreted counter-rotating gas loses angular momentum and flows inward, which lowers its measured rotation. That mechanism predicts the effect correlates with disturbance, which it does at 1 R_e but not at 1.5 R_e.
- Any follow-up on these data is **exploratory**, because the data have now been seen. A confirmatory test needs an independent sample (SAMI, CALIFA) or a newly registered discriminating statistic, such as the radial-flow harmonic s₁.

### 2026-10-02: RESULT ABOVE IS INVALID (position-angle convention bug)

The exploratory inflow analysis gave a median gas |V_R| of about 64 km/s, around 0.45 V_c, in *both* groups. That is physically implausible and exposed a bug.
- `pafit` returns PAs measured counter-clockwise from +y (that is, towards −x). The pipeline used them as if measured clockwise (towards +x).
- **Effect on ring fits:** they were done along the mirror image of the kinematic axis, misaligned by 2·PA (mod 180°). This biased every V_rot and leaked rotation into the radial term.
- **Effect on the receding-side flip:** it was correct only for axes within 45° of north–south. Galaxies near that boundary may have been misclassified between the counter-rotating and co-rotating groups.
- The synthetic validation did not catch it, because it supplied PAs directly and never used pafit's output.

**Commitment, made before re-running.**
1. The design is unchanged: amendment (g), with the same cuts, statistic, thresholds and verdict rules.
2. The convention is fixed: θ_pipeline = −θ_pafit.
3. A unit test now checks that the selection stage recovers a known receding PA from a synthetic velocity field.
4. Selection (ΔPA) and measurement are both re-run.
5. **The corrected run is the result of record, whatever it shows.**

The invalid outputs are archived in `results/invalid_pa_bug/` and will not be used. The earlier power checks also used the buggy geometry. They informed the choice of design but are not re-used as evidence.

### 2026-10-02: MaNGA G1 (design g), CORRECTED run. RESULT OF RECORD

PA convention fixed, with selection and measurement both re-run (see the bug note above). Log: `results/g1_manga_test_CORRECTED.log`.

**Sample.** After PA-quality and merger cuts: 112 counter-rotators (CR) and 4,166 co-rotators. After V_* > 40 km/s: 90 CR. Matched: 90 CR and 270 CO; the low-disturbance half has 46 CR.
- The old geometry inflated the CR count to 305 through random receding-side flips in weakly rotating galaxies.

**Primary test at 1 R_e.** S_pred = 42.2 km/s; the falsification threshold is 0.5 S_pred = 21.1 km/s.

| | k − 1 | k | k + 1 |
|---|---|---|---|
| S (km/s) | +7.9 ± 10.7 | +6.6 ± 8.6 | +14.2 ± 6.2 |
| S + 2σ | 29.3 | 23.8 | 26.6 |
| Low-disturbance half | +6.7 ± 12.9 | +8.3 ± 8.8 | +15.2 ± 7.4 |

**Registered verdict: INCONCLUSIVE.**
- Not SUPPORTED.
- Not FALSIFIED at the registered 50% bar: S + 2σ exceeds 21.1 at every k, by 2.7 km/s at the central k.

**Descriptive (not a registered verdict).**
- The measured slowdown is 16–34% of the prediction.
- The full predicted amplitude (S = S_pred) lies 3.2σ, 4.1σ and 4.5σ above the measurement at k − 1, k and k + 1.
- The 2σ upper bounds on a stellar-aligned flow are 69%, 56% and 63% of the RT-implied amplitude.

**Secondary results.**
- At 1.5 R_e: S = +16.6 ± 7.0, +19.4 ± 7.7 and +22.3 ± 7.3. The low-disturbance half gives +16 to +26.
- At 0.5 R_e: S ≈ 0.
- E is dominated by outliers and is uninformative.
- The co-rotators' median D at the central k is +11.3 km/s, a residual drift-calibration offset that cancels in S.

**Exploratory inflow analysis (corrected).** Median gas |V_R| is 6.7 km/s for CR and 5.8 for CO (about 4% of V_c), now physically plausible. D shows no correlation with |V_R|/V_c (|ρ| ≤ 0.18, p ≥ 0.09). For galaxies with |V_R| < 6 km/s, S = +7.3 ± 9.7 at the central k.

### 2026-10-02 (h): S2c details (registered before any black-hole mass is cross-matched)

**Context.** RT is no longer treated as a physical model. R_t remains an empirical descriptor: the radius where the rotation-curve excess saturates. S2c tests H_seed on that descriptor, using the existing Rational Taper fits.

**Black-hole masses.**
- Only dynamical measurements count: maser, stellar dynamics, gas dynamics, or reverberation mapping.
- Upper limits are excluded.
- Sources, in order of preference when a galaxy appears in more than one:
  1. van den Bosch 2016 (ApJ 831, 134) compilation;
  2. Kormendy & Ho 2013 (ARA&A 51, 511);
  3. Davis et al. 2017/2018 (spiral-galaxy M_BH).
- Matching is by name, including standard aliases (NGC/UGC/IC). Distances are rescaled to SPARC's: M_BH ∝ D.

**Sample.**
- Primary: the clean sample (Step 0).
- If the primary has fewer than 15 matches, the test is reported as EXPLORATORY. The taper-preferred sample (ΔBIC > 2, fit OK) is then also reported, labelled exploratory.

**Statistic.**
- The partial Spearman ρ(R_t/R_d, M_BH/M_bar | Π₂). Rank residuals are taken from linear fits on rank(Π₂).
- p-value from 10,000 permutations of the residuals.
- H_seed is supported only if p < 0.01, ρ > 0, and N ≥ 15.

**Also reported.**
- The same partial correlation with R_t alone, and with ω.
- The raw correlation with M_BH, as a check against the B/T proxy.

### 2026-10-02: S2c result (design h)

**Matching.**
- The SPARC × van den Bosch 2016 cross-match gives 12 galaxies.
- 6 are excluded as effective upper limits. Operational definition: lower uncertainty ≥ 1 dex, i.e. consistent with zero. Excluded: NGC 289, 300, 2903, 4088, 6503 and 7793.
- Davis+2017 adds no new SPARC galaxy.
- Remaining dynamical masses: NGC 3953, 3992, 4051, 5005, 5055 and 7331.

**Verdict: EXPLORATORY.** N_clean = 2, below the registered 15. No test of H_seed is possible.

**Exploratory numbers.**
- Partial ρ(R_t/R_d, M_BH/M_bar | Π₂): +1.00 (N = 4, taper-preferred); +0.97 (N = 6, all fits, permutation p ≈ 0.002).
- The sign matches H_seed.
- 4 of the 6 have unidentifiable R_t (they fail the clean-sample identifiability cut), so these numbers are not evidence.

### 2026-10-02: POST-HOC robustness check (not registered). RAR-normalised prediction

**Method** (`analysis/s_pred_rar.py`). The predicted S is recomputed by normalising the flow to the RAR instead of the registered SPARC flow-profile fits:
- u = V_RAR − V_bar, with g† = 1.2×10⁻¹⁰ m/s²;
- R·ζ = d(Ru)/dR;
- 123 SPARC galaxies (Q < 3, i ≥ 30°), matched as in the registered S_pred.

**Result.** S_pred,RAR = 33.0 km/s, against 42.2 km/s registered.

| | k − 1 | k | k + 1 |
|---|---|---|---|
| S_pred,RAR excluded at | 2.4σ | 3.1σ | 3.1σ |
| 2σ upper bound, as % of S_pred,RAR | 88% | 72% | 80% |

**Status.** This does not alter the registered verdict (INCONCLUSIVE). It is reported as a robustness check: the RAR normalisation is more conservative and gives weaker exclusion.

### 2026-10-02: POST-HOC robustness grid (not pre-specified)

**What was varied** (`analysis/g1_manga_robustness.py`; output `results/g1_manga_robustness.csv`):
- the speed cut V_* > 40, 60, 80 or 100 km/s;
- the radius, 1 R_e or 1.5 R_e;
- the drift factor, k − 1, k or k + 1;
- the amplitude normalisation, flow-profile template or RAR.

Matching and statistics are as in design (g).

**Finding.** S rises with the speed cut and with radius. At the central k:

| Speed cut | N_CR | S at 1 R_e (km/s) | S at 1.5 R_e (km/s) |
|---|---|---|---|
| > 40 | 90 | 6.6 ± 8.6 | 19.4 ± 7.7 |
| > 60 | 69 | 14.7 ± 9.5 | 27.1 ± 7.5 |
| > 80 | 58 | 13.7 ± 11.0 | 28.7 ± 8.1 |
| > 100 | 46 | 23.8 ± 13.5 | 32.0 ± 12.7 |

- The RAR-normalised prediction is 33–34 km/s. Of the four cuts, it is excluded at more than 3σ only at V_* > 40 and 1 R_e.
- At 1.5 R_e, or for V_* ≥ 60, the measurement is consistent with 50–100% of the RAR-normalised amplitude.

**Consequence.** The "disfavoured" statement depends on the primary analysis choice. The robust description is a positive asymmetry with the predicted sign, of 7–32 km/s, increasing with V_* and R. It is not distinguishable from gas disequilibrium with these data.

The pre-specified verdict (INCONCLUSIVE) stands. The subsamples are nested and not independent.

**Recommended confirmatory test.** An independent sample (SAMI or CALIFA), pre-registered on OSF before any data are downloaded, restricted to fast rotators (V_* > 80 km/s), with 1.5 R_e as the primary radius.

### 2026-10-02: EXPLORATORY discriminant analysis (informs the confirmatory pre-registration)

**Method** (`analysis/g1_manga_discriminant_exploratory.py`). For each counter-rotator, slowdown s_i = median(D_CO) − D_CR,i at the central k, with V_* > 40 km/s.

**Correlations.**
- With the vortex-predicted R·ζ_i (RAR): ρ = −0.10 (1 R_e) and −0.18 (1.5 R_e). The slowdown does not track the predicted R·ζ.
- With gas kinematic twist: ρ = +0.32 (p = 0.002, 1 R_e).
- With σ_gas/V_c: ρ = +0.22 (p = 0.04, 1.5 R_e).
- With radial flow: none.

**Radial growth.** The slowdown is +5.3 km/s at 1 R_e and +17.3 at 1.5 R_e, a ratio of about 3.3. The vortex-predicted ratio is 1.28.

**Twist split.** Median gas twist among the counter-rotators is 10.5°.

| | 1 R_e (km/s) | 1.5 R_e (km/s) |
|---|---|---|
| Low-twist (N = 44) | −6.8 ± 9.3 | +8.8 ± 6.7 |
| High-twist (N = 44 / 39) | +30.7 ± 12.4 | +30.0 ± 9.1 |
| Vortex prediction, both halves | ~29–31 | ~38–40 |

**Reading (exploratory, post hoc).**
- The asymmetry is carried by twisted, i.e. warped, unsettled or mis-deprojected, counter-rotating gas.
- Settled gas shows none. That is about 4σ below the vortex prediction at both radii.
- This favours a conventional explanation (disequilibrium, warps, or inclination mismatch in the deprojection) over a velocity-dependent field.

**Proposed confirmatory hypotheses for SAMI/CALIFA**, to be registered on OSF:
- **Primary:** the slowdown in low-twist, fast-rotating counter-rotators at 1.5 R_e.
  - The vortex reading predicts about R·ζ_pred.
  - The conventional reading predicts about 0.
- **Secondary:** the slowdown–twist correlation is positive.

### 2026-10-02: C1 confirmatory test registered on OSF

- **Registration:** https://osf.io/3dr29 (OSF Preregistration template). The text is in `preregistration/OSF_C1_confirmatory.md`.
- **Design.** The proposal in the entry above was refined before registration:
  - **Primary:** the Theil-Sen zero-twist intercept of the slowdown at 1.5 R_e, judged by a Bayes factor between H_V (a = P, the RAR-normalised prediction) and H_C (a = 0). The fraction f = a/P is reported as the estimand.
  - **Secondary:** analyses 1–9 as listed in the registration.
- **Frozen code:** commit f3cc6de9aa2296abdc0736f07143acb9deaea0a9 (tag `c1-prereg`), archived at https://doi.org/10.5281/zenodo.23112005. After the tag, only the text of the OSF form changed.
- **Status of the data:** no map of any galaxy eligible for C1 has been accessed. The SAMI and CALIFA sample maps are downloaded only after the registration is approved.
- **Allowed changes from here on:** only mechanical loader fixes found in the dry run, each logged here with a date before the full run.

### 2026-10-02: C1 registration approved; data access begins

- OSF registration https://osf.io/3dr29 was approved by the author.
- Only after that were the C1 sample maps requested: the SAMI bulk download (Data Central) and the CALIFA download (`analysis/c1_fetch_califa.py`, which only downloads, with no measurement and no output beyond progress counts).

### 2026-10-02: C1 loader fix (mechanical, before the select stage)

- **Problem.** The SAMI bulk download contains repeat-observation cubes (B, C, D) for 180 galaxies. The frozen reader required a unique file match per product, so it would have treated those galaxies as missing.
- **Fix.** `src/sami.py` now selects the ISBEST cube from CubeObs, as the registration already specifies ("ISBEST = 1 selects one cube per galaxy"), and indexes the map directory once. A unit test covers the change.
- **Result.** All 2141 SAMI parent galaxies have complete map sets. 34 of them use a best cube other than A.
- **What did not change:** cuts, statistics and decision rules. No measurement had been run.

### 2026-10-02: C1 select and dry run (counts only)

**Position-angle stage.**
- 2391 parent galaxies, with no failures. SAMI is 2141 of these and CALIFA 250.
- 728 galaxies pass the PA-quality cuts.
- Counter-rotator candidates: 8 (SAMI 7, CALIFA 1).
- Co-rotating candidates: 678.

**Why the yield is low.** The MaNGA-derived rule of at least 50 valid spaxels after 2× decimation fails for gas in 994 of 2141 SAMI galaxies (SAMI's field of view is 15″) and for stars in 843. Among the 15 raw ΔPA > 150° candidates, PA errors above 20° remove about half. The registered yield estimate of 15–20 was based on published counter-rotator counts and did not apply these rules to SAMI's field of view. That was an error in the sample-size forecast, not a change of plan.

**Dry run.**
- All 8 CR pass V_* > 40 and are matched 3:1, giving 24 controls.
- Only 5 CR (all SAMI) have valid 1.5 R_e measurements and a twist.
- Under the registered stopping rule (fewer than 10), inference is not carried out, and the run reports descriptively.
- The dry run required no loader fixes.

**Noted before the full run.** The registration lists the secondary analyses as reported "regardless of outcome", but the frozen `c1_test.py run` computes them only when N ≥ 10. After the registered run, the secondary quantities will be computed with the same frozen functions and reported as descriptive. No inference will be drawn from them, since N < 10.

### 2026-10-02: C1 registered run (executed once). RESULT OF RECORD

**Registered verdict: UNDERPOWERED (N = 5 < 10), descriptive only.** All 5 primary counter-rotators are in SAMI. The single CALIFA counter-rotator, NGC 1167, has no valid 1.5 R_e ring.

**Descriptive quantities** (`results/c1.json`, `results/c1_descriptive.json`):
- Prediction: P = 44.9 km/s.
- Per-CR slowdown s at 1.5 R_e, central k, with the twist in brackets:
  - 230421: +3.3 (2°)
  - 300787: +3.4 (1°)
  - 278840: +6.6 (0°)
  - 321059: −5.4 (47°)
  - 184648: +172.6 (6°)
- Theil-Sen intercept: a = +3.7 ± 38.1 km/s, so f = 0.08. The likelihood ratio is 0.56 (it would read "inconclusive" even if N were larger). The same holds at k ± 1.
- Four of the five sit within 7 km/s of zero. One outlier (184648) inflates σ_a.
- Secondary 1: S at 1.5 R_e = +27.8 ± 13.7 km/s (central k; 6 CR, 17 controls). This is close to the MaNGA value. The statistic mixes settled and unsettled gas.
- Secondary 9 (balance): the matching covariates have standardised mean differences below 0.25. The CR gas is hotter at 1.5 R_e (σ_g 57 against 43 km/s, SMD 0.64). V_c,* is somewhat lower in the CR (SMD −0.33).

**Secondary 7.** All 5 SAMI primary CR are also counter-rotators (ΔPA > 150°) in the SAMI DR3 catalogue PAs. The pipeline's classification is confirmed. N < 10, so the primary statistic was not computed.

**Secondary 5 (pooled MaNGA exploratory + C1; NOT confirmatory).**
- N = 85.
- a = +7.6 ± 6.3 km/s against P = 40.6, so f = 0.19 with a 95% upper bound of 0.55.
- log10 likelihood ratio = −5.6, favouring H_C.
- This is dominated by the MaNGA data, which generated the hypothesis, so it is reported as exploratory.

**Secondary 8:** running.

### 2026-10-02: C1 secondary 8, and a correction to how secondary 5 is read

**Secondary 8 (MaNGA, with the registered per-CR control definition of s).**
- With the merger cut: a = +15.6 ± 10.8 km/s, P = 40.3, f = 0.39 (95% CI −0.02 to 0.96; one-sided upper bound 0.90). The likelihood ratio is 0.21, labelled "moderate evidence for H_C". N = 83.
- Without the merger cut: a = +18.6 ± 9.7, f = 0.46 (upper bound 0.83). The likelihood ratio is 0.51, labelled "inconclusive". N = 83.
- So the merger cut has a small effect.

**Correction (found while reading secondary 8).** The exploratory MaNGA file `g1_manga_discriminant_cr.csv` defines s_i = median(D over *all* matched controls) − D_i. The C1 registration defines s_i using *each CR's own* 3 controls. The frozen secondary 5 pools MaNGA s values of the first kind with C1 values of the second, so its result (f upper bound 0.55) mixes definitions. It should not be quoted as a constraint.

The consistent MaNGA figure is secondary 8 with the merger cut: f upper bound 0.90. The C1 power analysis (`c1_power.py`) was calibrated on the global-median definition, so its residual scatter, and hence its power, were probably optimistic.

### 2026-10-03: POST HOC comparison with the strong-gravitomagnetic GR model (not registered)

**Model.** Astesiano, Ruggiero & Re (2026, arXiv:2606.18868), Eq. 17. Their Eq. 4 is a geodesic equation, so counter-rotating gas follows the counter-rotating branch. The derivation is `derivations/04_strong_gravitomagnetic_branches.py`. The predicted asymmetry is ΔV = (V_obs² − V_bar²)/V_obs for their solution ψ = C0 r, and ΔV = 2(V_obs − V_bar) for rigid dragging. It is evaluated with the RAR on matched SPARC galaxies, at each counter-rotator's own radius.

**Statistic.** The C1 primary statistic (zero-twist intercept, with per-CR controls), applied to MaNGA with P replaced by the model prediction. Code: `analysis/agm_comparison.py`; output in `results/agm_comparison.*`.

| Radius | Model ΔV | a (central k) | f, 95% upper | log10 LR (k−1 / k / k+1) |
|---|---|---|---|---|
| 1.0 R_e | 32 | −6.9 ± 8.8 | 0.21 | −4.1 / −4.1 / −1.6 |
| 1.5 R_e | 40.5 | +11.6 ± 10.6 | 0.76 | −1.4 / −1.4 / +0.8 |

- At 1.0 R_e the model's counter-rotating branch is strongly disfavoured at all three drift factors.
- At 1.5 R_e it is disfavoured at k−1 and k, but moderately favoured at k+1. The drift correction is the controlling systematic there.
- The rigid-dragging variant gives similar or stronger results.
- SAMI (descriptive, N = 5): model 45 km/s; observed s = +3.3, +3.4, +6.6, +172.6, −5.4.

**Caveats.**
- The authors apply ψ = C0 r beyond the luminous disc. These radii are inside it, so this tests an extension of the model.
- The MaNGA data are exploratory (already seen).
- This rebuild of the MaNGA sample gives a slightly different 1.5 R_e intercept (11.6 against 15.6 in secondary 8, with N = 80 against 83) because the measurement cache and inclusion rules differ.

**Outer-disc forecast** for a typical counter-rotator (V_c ≈ 172 km/s, R_e ≈ 1.7 kpc):
- Predicted ΔV: 45, 60 and 73 km/s at 2, 3 and 4 R_e.
- A decisive (likelihood ratio > 10, 80% power) zero-twist test needs about 40–80 counter-rotators at 3 R_e, or 60–120 at 2 R_e, depending on the measurement scatter.


### 2026-10-04: CORRECTION to the flowing-space prediction; full number audit of the manuscript

**Error found.** The "flowing space" prediction used since the 2026-10-02 RAR
normalisation (`analysis/s_pred_rar.py`, `pred_rz`, the C1 registration's P and the
"FS" rows of the draft Table 1) was ΔV = d(Ru)/dR with u = V_RAR − V_bar. That assumed
the additive relation V_obs = V_bar + u, which holds only for solid-body flow. For a
general flow the prograde root of Formulation A is v_pro = V_bar + u − SR/2 + O(S²R²/V_b)
(`derivations/01_force_laws.py` [A6], new check), so the first-order-consistent
prediction from an observed rotation curve is ΔV = 2(V_obs − V_bar) for any flow profile,
the same as Formulation B and as rigid dragging (ψ ∝ r²). The old formula gives
32.1 / 40.3 km/s at 1.0 / 1.5 R_e for the MaNGA counter-rotators; the corrected
value is 36.8 / 46.8 km/s. Solving Formulation A exactly for u(R)
(`analysis/paper_checks.py` `fs_shear`, an ODE) differs from 2(V_RAR − V_bar) by at most
1.9 km/s per counter-rotator.

**Consequences.**
- The paper's Table 1 now has two predictions per radius: ψ ∝ r (SGM1; unchanged) and
  flowing space / vector potential / rigid dragging (FS = the former SGM2 row). The former
  FS row (32.1 / 40.3 km/s) is withdrawn. The headline numbers (a, f < 0.21, the
  likelihood-ratio ranges, the 32–37 and 40–47 km/s prediction ranges) are unchanged.
- The robustness checks (merger cut, matching seeds, pooled-control definition) now quote f
  and the likelihood ratio against the smaller ψ ∝ r prediction (P = 32.0 instead of 32.1
  km/s at 1 R_e); the changes are in the second decimal.
- The C1 registered P (44.9 km/s for the five SAMI counter-rotators) was computed with the
  superseded formula. It is retained for the registered analysis (which was underpowered
  in any case) and happens to equal the ψ ∝ r value (44.8); the corrected flowing-space
  value is 51 km/s. The paper states this.
- `results/paper_numbers.json` keeps the superseded value under `pred_FS_superseded` and
  `SAMI.P_registered_FS_superseded`; `analysis/s_pred_rar.py` and `pred_rz` are unchanged
  so that the registered C1 numbers remain reproducible.

**Other corrections from the same audit** (manuscript only; all numbers re-derived from
`results/` and found to match): the sign of the vector-potential convention in Eq. (2)
((∇×A)_z = −2Ω, as in `derivations/01`); the abstract's likelihood-ratio statement now
restricted to the central and lower drift factors; `paper_checks.py` counted duplicated
rows of `manga_pa.csv` (14464) for the "all 7582 measured" claim, fixed to count galaxies;
several roundings tightened. A new symbolic proof replaces the grid check for
"ΔV ≥ V_obs − V_bar for ψ ∝ rⁿ" (`derivations/04`; it holds for all n ≥ 1).
