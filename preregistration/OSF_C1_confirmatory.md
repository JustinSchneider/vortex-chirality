# OSF Preregistration: confirmatory chirality test C1 (SAMI DR3 + CALIFA DR3)

DRAFT, 2026-10-02. Sections follow the fields of the OSF "OSF Preregistration"
template, in order. Items marked [[...]] await the data-product documentation
check and must be resolved before submission.

---

## Title

Do gas disks counter-rotating with respect to their stars rotate more slowly? A
pre-registered confirmatory test with SAMI and CALIFA integral-field kinematics

## Description

Some proposed alternatives to dark matter add a velocity-dependent force to
galaxy dynamics. Two examples are vorticity of space (a "flowing space" or
gravitomagnetic-like field) and a Coriolis-type term. Such a force is odd under
time reversal, so it treats orbits that move with and against the field
differently. At a given radius, stars and gas co-rotating with the field gain
circular speed. Material rotating the other way loses it. The predicted
asymmetry is ΔV = Rζ, where ζ is the field's vorticity. Dark matter, MOND and
modified inertia are time-reversal symmetric and predict ΔV = 0 after correcting
for pressure support (asymmetric drift).

Galaxies whose ionised gas counter-rotates with respect to their stars are a
natural laboratory. The field is set by the dominant (stellar or baryonic)
rotation. Under a velocity-dependent term, the counter-rotating gas should show
a lower drift-corrected circular speed than the stars, relative to matched
co-rotating controls.

An exploratory analysis of MaNGA DR17 (90 counter-rotators, 270 matched
controls; full log at [[GitHub/Zenodo URL]]) gave four findings:
- The counter-rotators show a mean "slowdown" that is smaller than the full
  vortex prediction.
- The slowdown is not robust to analysis choices.
- The slowdown correlates with the kinematic twist of the gas (Spearman
  ρ = +0.32, p = 0.002), not with the predicted vortex amplitude.
- When the slowdown is extrapolated to zero gas twist at 1.5 R_e, it is
  +11.6 ± 7.4 km/s, against a vortex prediction of about 40 km/s.

Those results were found after looking at the data, so they generate the
hypotheses tested here; they do not confirm them. This registration fixes, in
advance, a test of the zero-twist slowdown on independent samples: SAMI DR3 and
CALIFA DR3, excluding any galaxy also in MaNGA.

## Hypotheses

Let s_i be the slowdown of counter-rotator i at 1.5 R_e: the median D of its
matched co-rotating controls minus its own D, where D is the drift-corrected
gas-minus-stars circular-speed difference (defined under "Indices"). Let t_i be
its gas kinematic twist in degrees. Fit s_i = a + b·t_i. The intercept a is the
slowdown of a counter-rotator whose gas has settled into a single plane.

- **H_V (velocity-dependent / vortex term carries the mass discrepancy):**
  a = P, where P is the mean over the counter-rotators of the predicted Rζ_i.
  It is computed by assuming that the field accounts for all of the
  radial-acceleration-relation (RAR) discrepancy (see "Statistical models").
  Directional: a > 0. Expected P ≈ 35–45 km/s, based on the MaNGA sample.
- **H_C (conventional: dark matter / MOND / modified inertia, with
  non-equilibrium gas):** a = 0, with b > 0, i.e. slowdown grows with gas twist
  because unsettled, warped or inflowing gas is not on circular orbits.

The primary estimand is the fraction f = a / P. This is the share of the full
vortex amplitude that the data allow. H_V is f = 1; H_C is f = 0.

## Study type

Observational study, analysing existing public data.

## Blinding

There are no human participants. Analysis blinding works as follows:
- The entire pipeline (selection → measurement → matching → statistics) is
  frozen at the commit given below before any SAMI or CALIFA kinematic product
  is downloaded.
- A "dry" mode will first be run to check mechanics. It prints only sample
  counts. It does not print D, s, a, b, f or any rotation amplitude.
- The survey loaders are written from the documentation alone. During the dry
  run, a loader may be fixed only for a mechanical error: a crash, a
  mis-read extension or unit, or a wrong mask convention. Every fix is
  committed and logged with its date in `preregistration/PREREG.md` before the
  full run. No change to cuts, statistics or decision rules is allowed.
- The full analysis is then run once.

## Study design

This is a case-control comparison, run in the same way on each survey:

1. **Measure kinematic position angles.** Every galaxy with both stellar and
   ionised-gas velocity maps gets stellar and gas PAs, measured with the
   repository's `receding_pa` (pafit, Krajnović et al. 2006, with the
   convention fix and unit test in commit 1ee3136).
2. **Form the groups.**
   - Counter-rotators (CR): ΔPA > 150°.
   - Co-rotating pool (CO): ΔPA < 30°.
   - Both groups also require PA errors ≤ 20° for both components.
3. **Match controls.** Each CR gets 3 CO controls, chosen by greedy
   nearest-neighbour matching without replacement on standardised covariates.
   Seed 20261002, the same algorithm as the MaNGA run.
4. **Measure rings.** All CR and controls are measured in elliptical rings at
   0.5, 0.75, 1.0, 1.25 and 1.5 R_e, each 0.25 R_e wide.

## Randomization

Not applicable. The order of greedy matching is a fixed-seed random
permutation (seed 20261002).

## Existing data

Registration prior to accessing the data.

## Explanation of existing data

SAMI DR3 and CALIFA DR3 are public. As of registration, the author has not
downloaded or viewed any SAMI or CALIFA stellar or gas velocity map, cube, or
catalogue value of rotation velocity or dispersion.

The author has read:
- the SAMI DR3 and CALIFA DR3 release papers;
- the SAMI misalignment study of Ristea et al. (2022), which reports counts of
  misaligned and counter-rotating galaxies but no comparison of gas and
  stellar circular speeds;
- the CALIFA kinematic-alignment studies of Barrera-Ballesteros et al. (2014,
  2015) and García-Lorenzo et al. (2015), which report position angles and
  alignment, but not gas-versus-stellar circular-speed differences;
- the SAMI DR3 and CALIFA data-product documentation, including table and
  column names (but no values).

**Incidental exposures, disclosed for completeness.** Neither one reveals
rotation or circular-speed amplitudes.
1. While checking the structure of Ristea et al. (2022) Table D1, the first
   four data rows were displayed. They show the CATID, stellar and gas PA,
   SFR, stellar mass, and misalignment cause.
2. The SAMI documentation's example-query page states the number of rows that
   a slow-rotator query returns.

The author has also analysed MaNGA DR17 with the same pipeline. Those
exploratory results motivate the hypotheses and the power analysis.

The pipeline measures every kinematic PA itself. The published PAs (the SAMI
catalogue PA_STELKIN and PA_GASKIN, Ristea et al. Table D1, and the
Barrera-Ballesteros tables) are not used for selection. They are used only in
secondary analysis 7.

Galaxies observed by MaNGA are excluded, so no counter-rotator can be in both
samples.

## Data collection procedures

**SAMI DR3** (Croom et al. 2021). Maps are retrieved through the Data
Central bulk-download service; catalogues come through its TAP service.
- **Stellar maps:** stellar velocity and dispersion. The two-moment pPXF fit,
  default (unbinned) 0.5″ spaxels.
  - Extensions: VEL, VEL_ERR, SN; and SIG, SIG_ERR.
  - Valid velocity spaxels: VEL_ERR < 30 km/s and SN > 3.
  - Valid dispersion spaxels: SIG_ERR < 0.1 SIG + 25, SN > 3, SIG > 35 km/s,
    and VEL_ERR < 30 km/s. These are the SAMI-recommended cuts (van de Sande
    et al. 2017; Croom et al. 2021).
- **Gas maps:** ionised-gas velocity and dispersion. The LZIFU "1-comp"
  (single-Gaussian) fit, which ties the velocity across all strong lines.
  - Extensions: primary plus V_ERR; and primary plus VDISP_ERR.
  - Valid gas spaxels: V_ERR ≤ 30 km/s and Hα S/N ≥ 5, using the 1-comp Hα
    flux and its error.
- **Catalogues:**
  - `sami_dr3.InputCatGAMADR3` and `InputCatClustersDR3` (whichever
    applies): z, stellar mass, BAD_CLASS.
  - `sami_dr3.MGEPhotomUnregDR3`: R_e = ReMGE, ε = epsMGE_Re, photometric
    PA = PAMGE.
  - `sami_dr3.VisualMorphologyDR3`: TYPE.
  - `sami_dr3.CubeObs`: ISBEST = 1 selects one cube per galaxy.
- **Not used:** the SAMI kinematic tables (`samiDR3Stelkin`,
  `samiDR3gaskinPA`) are not queried before the full run.

**CALIFA** (Sánchez et al. 2016).
- **Stellar kinematics:** the V1200 products of Falcón-Barroso et al.
  (2017), about 300 galaxies: `<name>.CALIFA.V1200.stekin.fits`.
  - These are binary tables of Voronoi-binned pPXF fits. Each spaxel (X, Y)
    takes the velocity Vp ± DVp and the dispersion Sp ± DSp of its bin.
  - Valid spaxels: QC flag good; DVp ≤ 30 km/s.
- **Gas kinematics:** the eCALIFA pyPipe3D V500 products (Sánchez et al.
  2023), ELINES extension.
  - Channel 0 is the Hα velocity.
  - Channel 1 is the Hα FWHM in Å, converted to an intrinsic dispersion:
    σ = c·√[(FWHM / 2.354 / λ_Hα)² − (2.6 Å / λ_Hα)²].
  - Hα flux and its error come from FLUX_ELINES.
  - Valid gas spaxels: Hα S/N ≥ 5; velocity error ≤ 30 km/s.
- **Catalogue quantities:** R_e, ellipticity, photometric PA, stellar mass,
  and Hubble type, all from the CALIFA DR3 / eCALIFA catalogues
  [[exact table, confirmed during loader writing]].
- **Sample:** only galaxies with both a stellar and a gas product enter.
- **Gas–stellar map registration:** the two maps are put on a common frame
  using the WCS and header reference pixel (XCEN, YCEN for the stellar
  product). If they cannot be registered to within 1″, the galaxy is
  excluded.

**Instrumental dispersion.** The pPXF stellar dispersions (both surveys) and
the LZIFU gas dispersions (SAMI; Zhou et al. 2017) are taken to be intrinsic,
i.e. instrument-corrected. Any residual error in this correction enters
counter-rotators and their matched controls alike. Because it appears only
through the k σ² term, it largely cancels in s.

**Overlap exclusion:**
- Any SAMI or CALIFA galaxy within 3″ of a MaNGA DR17 target is removed. The
  published overlaps are about 74 SAMI–MaNGA galaxies (Fraser-McKelvie et
  al. 2021) and a few tens of CALIFA–MaNGA galaxies.
- Galaxies in both SAMI and CALIFA are kept once, from CALIFA, because its
  wider field covers 1.5 R_e more often.

## Sample size

Every galaxy that passes the cuts is analysed; there is no target N. Based on
published counter-rotator counts and coverage, the expected yield is 15–20
counter-rotators with valid measurements at 1.5 R_e.

## Sample size rationale

The power analysis is `analysis/c1_power.py`; its output is
`results/c1_power.json`. It resamples the MaNGA twist distribution and
Theil-Sen residuals (robust SD 46 km/s). It holds b at its MaNGA value under
both hypotheses and uses P = 40 km/s.

The table gives the probability that the Bayes factor favours the true
hypothesis by more than 10, and, after the slash, the probability that it
favours the wrong one by more than 10:

| N_CR | H_C true | H_V true |
|---|---|---|
| 10 | 0.21 / 0.005 | 0.12 / 0.018 |
| 15 | 0.33 / 0.002 | 0.28 / 0.018 |
| 20 | 0.50 / 0.002 | 0.45 / 0.018 |
| 25 | 0.60 / 0.002 | 0.56 / 0.027 |
| 30 | 0.72 / 0.001 | 0.66 / 0.022 |

At the expected N, a decisive result is about as likely as not. A wrong
decisive result has a probability of 3% or less. This limitation is accepted
and stated in advance. An inconclusive outcome will be reported as such,
together with the interval on f.

## Stopping rule

The analysis is run once on the full released samples. If fewer than 10 CR have
valid 1.5 R_e measurements, inference is not carried out. The primary
quantities are then reported descriptively, and the result is declared
underpowered.

## Manipulated variables

None.

## Measured variables

Per galaxy, all derived from survey products by the frozen pipeline:
- stellar and gas kinematic PA, and ΔPA;
- inclination from the photometric axis ratio, with intrinsic thickness
  q0 = 0.2, and i ≥ 30° required;
- R_e (arcsec, and kpc from the redshift distance);
- stellar mass;
- morphology on a coarse common scale from 0 to 3:
  - SAMI: TYPE as published (0 E, 0.5 E/S0, 1 S0, 1.5, 2 early spiral, 2.5,
    3 late spiral/irregular). TYPE 5 or −9 counts as missing.
  - CALIFA: Hubble type mapped as E = 0, S0 = 1, Sa–Sb = 2, Sbc and later
    = 3.
  - Morphology is used only for matching, which happens within a survey.
- σ_* within 1 R_e, as the median of instrument-corrected stellar dispersion in
  spaxels inside the 1 R_e ellipse;
- in each ring, by harmonic fit (V_sys + c1 cos φ + s1 sin φ, requiring at least
  3 azimuthal quadrants):
  - stellar and gas V_rot and |V_R|;
  - the median stellar and gas dispersion;
  - the gas residual RMS;
- gas twist t: the PA difference between gas inside and outside 1 R_e, each
  needing at least 30 valid spaxels;
- the drift factor k = R/h_R + R/h_σ² − ½, with h_R = R_e/1.678, h_σ² from a
  log-linear fit to the ring dispersions, and k clipped to [0, 4].

**Spaxel validity:**
- survey quality mask clean;
- velocity error ≤ 30 km/s;
- gas: Hα S/N ≥ 5;
- dispersions: instrument-corrected; gas below instrumental resolution set
  to 0.

**Ring validity:** valid spaxels must cover at least 5 arcsec² (20 SAMI
spaxels of 0.5″, or 5 CALIFA spaxels of 1″) and at least 3 azimuthal
quadrants. This matches the MaNGA rule of 20 spaxels of 0.5″. For CALIFA
stellar maps, a Voronoi bin counts once per spaxel it covers.

## Indices

- D = (V_g² + k_g σ_g² − V_*² − k_* σ_*²)/(V_g + V_*), in km/s, at each ring.
- s_i = median(D over the 3 controls of i) − D_i.
- P = mean over the CR of Rζ_i. Rζ_i is evaluated at 1.5 R_e,i (kpc) for each
  SPARC galaxy (Q < 3, i ≥ 30°) with |V_flat − V_c,*,i| < 30 km/s, and then
  averaged:
  - Rζ = d(R u)/dR;
  - u = V_RAR − V_bar;
  - the RAR uses g† = 1.2×10⁻¹⁰ m s⁻²;
  - V_c,* = √(V_*² + k_* σ_*²) at 1 R_e.
  - The code is `analysis/s_pred_rar.py` (`rar_profiles`) together with
    `pred_rz` in `analysis/g1_manga_discriminant_exploratory.py`, unchanged.
- f = a / P.

## Statistical models

**Primary.** The model is a Theil-Sen regression of s_i on t_i over all CR with
valid 1.5 R_e measurements, with the central drift factor k:
- The intercept is a = median(s − b t), where b is the median of the pairwise
  slopes.
- σ_a is the SD of a over 10,000 bootstrap resamples of the CR, with seed
  20261002. Controls stay attached to their CR in each resample.
- The Bayes factor of H_V against H_C is

  BF = exp[(a² − (a − P)²) / (2 σ_a²)].

- The point estimate and 95% bootstrap percentile interval of f = a / P are also
  reported.

**Matching covariates** (standardised):
- log M_*;
- morphology (0–3 scale);
- log σ_*(1 R_e);
- inclination;
- survey (exact match: SAMI CR are matched to SAMI controls, and CALIFA CR to
  CALIFA controls).

**Secondary analyses.** These are pre-specified and reported regardless of
outcome, but none of them decides the verdict:
1. As in the MaNGA run, S = median(D_CO) − median(D_CR) at 1.0 and 1.5 R_e,
   with bootstrap SD, at k−1, k and k+1.
2. b > 0, one-sided. This tests the conventional prediction that slowdown rises
   with twist. The test is a bootstrap 95% lower bound on b.
3. The primary BF and f at k−1 and k+1.
4. The primary analysis for each survey separately.
5. A pooled MaNGA + SAMI + CALIFA analysis. It is labelled as combining
   exploratory and confirmatory data.
6. The Spearman correlation of s_i with predicted Rζ_i. H_V predicts a
   positive value.
7. The agreement of the pipeline's CR classification with published PAs: the
   SAMI catalogue PA_STELKIN/PA_GASKIN, Ristea et al. (2022) Table D1, and
   Barrera-Ballesteros et al. (2014, 2015). The primary analysis is then
   repeated on the CR on which both agree.
8. The MaNGA exploratory primary statistic (intercept a, BF, f), recomputed
   without the Galaxy Zoo merger cut. This shows the size of the effect of
   omitting the merger cut.

## Transformations

None beyond the indices above. Twist is in degrees, from 0 to 180.

## Inference criteria

The following labels apply to the primary analysis at central k:

| Bayes factor | Label |
|---|---|
| BF ≥ 10 | Strong evidence for H_V |
| 3 ≤ BF < 10 | Moderate evidence for H_V |
| 1/10 < BF ≤ 1/3 | Moderate evidence for H_C |
| BF ≤ 1/10 | Strong evidence for H_C |
| Between 1/3 and 3 | Inconclusive |

- If the label at k−1 or k+1 is weaker than at central k, the weaker label is
  reported as the headline.
- Regardless of the label, the 95% upper bound on f is reported as a
  constraint on any velocity-dependent contribution.

## Data exclusion

- Galaxies are excluded if they:
  - are in MaNGA DR17;
  - lack either velocity map;
  - have a PA error above 20°;
  - have i < 30°;
  - have V_*(1 R_e) ≤ 40 km/s;
  - fail the survey's catalogue-level quality flags (SAMI: a BAD_CLASS
    value marking a bad target [[values from InputCat docs]], or not
    ISBEST; CALIFA: the
    survey's own exclusion flags).
- Neither survey has a merger or disturbance flag equivalent to Galaxy Zoo's,
  so no merger cut is applied. This is a deliberate deviation from the
  MaNGA run. Secondary analysis 8 measures its effect.
- CR without a valid 1.5 R_e ring (for both stars and gas) or without a valid
  twist are excluded from the primary analysis. They are counted in the
  report.
- There is no outlier removal; Theil-Sen is used for its robustness.

## Missing data

Missing catalogue covariates exclude a galaxy from matching. Missing rings
exclude a galaxy only from analyses at that radius.

## Exploratory analysis

Anything not listed above will be labelled exploratory. That includes
alternative radii, alternative velocity cuts, alternative matching, and models
in which a depends on stellar mass.

## Other

**Code.**
- Repository: [[GitHub URL]].
- Frozen commit: [[hash]], tagged `c1-prereg`.
- Archive: [[Zenodo DOI]].
- Any deviation from this plan will be listed, dated and justified in
  `preregistration/PREREG.md` and in the paper.

**Disclosure.** The analysis code and this plan were developed with an AI coding
assistant (Claude, Anthropic). All decisions are the author's.
