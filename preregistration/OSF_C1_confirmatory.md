# OSF Preregistration: confirmatory chirality test C1 (SAMI DR3 + CALIFA DR3)

DRAFT, 2026-10-02. Sections follow the fields of the OSF "OSF Preregistration"
template, in order. Items marked [[...]] must be filled in before submission.

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
controls; full log at https://github.com/JustinSchneider/vortex-chirality) gave four findings:
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
3. **Match controls.** Each CR gets 3 CO controls (fewer only under the
   fallback rule in "Statistical models"), chosen by greedy nearest-neighbour
   matching without replacement on standardised covariates.
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
3. The eCALIFA catalogue files contain integrated kinematic columns (Vmax,
   vel_sigma_Re, Lambda_Re, and others). The files were downloaded and their
   column names listed, but those columns were never read, printed or used.
4. For the format check, maps of galaxies that are also in MaNGA (and
   therefore excluded from C1) were downloaded and inspected: CALIFA
   UGC 08107 and SAMI 517164.

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
  - Files `<CATID>_A_stellar-velocity_default_two-moment.fits` and
    `..._stellar-velocity-dispersion_default_two-moment.fits`. Extensions:
    VEL, VEL_ERR, SNR; and SIG, SIG_ERR, SNR.
  - Valid velocity spaxels: VEL_ERR < 30 km/s and SNR > 3.
  - Valid dispersion spaxels: SIG_ERR < 0.1 SIG + 25, SNR > 3, SIG > 35 km/s,
    and VEL_ERR < 30 km/s. These are the SAMI-recommended cuts (van de Sande
    et al. 2017; Croom et al. 2021).
- **Gas maps:** ionised-gas velocity and dispersion. The LZIFU "1-comp"
  (single-Gaussian) fit, which ties the velocity across all strong lines.
  - Files `..._gas-velocity_default_1-comp.fits` (primary plus V_ERR),
    `..._gas-vdisp_default_1-comp.fits` (primary plus VDISP_ERR) and
    `..._Halpha_default_1-comp.fits` (plane 0, the total flux, plus
    HALPHA_ERR). All are 50×50 maps with 0.5″ spaxels, centred on the
    galaxy (CRPIX 25.5).
  - Valid gas spaxels: V_ERR ≤ 30 km/s and Hα S/N ≥ 5, using the 1-comp Hα
    flux and its error.
- **Catalogues** (TAP; non-kinematic columns only; copies in `data/sami/`):
  - `sami_dr3.InputCatGAMADR3` and `InputCatClustersDR3`, whichever applies:
    - position;
    - redshift: z_tonry (flow-corrected) where available, otherwise z_spec;
    - R_e = r_e (Sérsic, major axis, as for MaNGA);
    - ellipticity `ellip`, with b/a = 1 − ellip;
    - stellar mass `Mstar`;
    - BAD_CLASS.
  - `sami_dr3.VisualMorphologyDR3`: TYPE.
  - `sami_dr3.CubeObs`: ISBEST = 1 selects one cube per galaxy; warning flags.
  - Angular scale: flat ΛCDM with H0 = 70 and Ωm = 0.3.
- **Quality.** A galaxy is kept only if:
  - BAD_CLASS ∈ {0, 5, 8} (the "good target" values; Croom et al. 2021,
    Sect. 2.4);
  - WARNSTAR, WARNSK2M, WARNSKER, WARNWCS, WARNEMFT and WARNMULT are all 0.
- **Map orientation.** Maps are taken as centred on the galaxy (the WARNWCS
  galaxies are excluded).
- **Locating files.** DR3 file names are not documented, so each product is
  located by keywords in the file name, and the match must be unique (code:
  `src/sami.FILE_KEYS`).
- **Not used:** the SAMI kinematic tables (`samiDR3Stelkin`,
  `samiDR3gaskinPA`) are not queried before the full run.

**CALIFA** (Sánchez et al. 2016; eCALIFA, Sánchez et al. 2023).
- **Stellar kinematics:** the V1200 products of Falcón-Barroso et al.
  (2017), 300 galaxies: `<name>.CALIFA.V1200.stekin.fits`.
  - These are binary tables of Voronoi-binned pPXF fits, with one row per 1″
    spaxel. X and Y are in arcsec, relative to the galaxy centre.
  - Each spaxel carries the velocity Vp ± DVp and the dispersion Sp ± DSp of
    its bin.
  - Valid velocity spaxels: DVp ≤ 30 km/s.
  - Valid dispersion spaxels: additionally DSp < 0.1 Sp + 25 and
    Sp > 35 km/s. This mirrors the SAMI rule.
  - The QC column is a continuous fit-residual measure, not a flag, so it is
    not used.
- **Gas kinematics:** the eCALIFA pyPipe3D V500 v2.3 cubes, 0.5″ spaxels.
  - FLUX_ELINES planes "vel Ha", "e_vel Ha", "flux Ha" and "e_flux Ha" give
    the velocity, its error, the flux and its error.
  - ELINES channel 1 is the Hα FWHM in Å, including instrumental width. It is
    converted to an intrinsic dispersion:
    σ = c·√[(FWHM / 2.354 / λ)² − (2.6 Å / λ)²], with λ = 6562.8 (1+z) Å.
    Values below the instrumental width are set to 0.
  - Valid gas spaxels:
    - Hα S/N ≥ 5;
    - 0 < velocity error ≤ 30 km/s;
    - outside the Gaia foreground-star mask.
- **Catalogue quantities** come from the eCALIFA tables
  (`galaxies_properties.fits`: position, z, Hubble type, QC_flag;
  `eCALIFA.pyPipe3D.fits`: Re_arc, ellip, log_Mass, DA). Two columns are not
  what their names suggest; both readings are verified on all 895 rows:
  - `ellip` is an eccentricity, since it tracks `ecc` with r = 0.998. So
    b/a = √(1 − ellip²).
  - `DA` is the angular scale in kpc/arcsec, since it equals Re_kpc / Re_arc.
  - Galaxies with QC_flag = 1 are excluded.
- **Sample:** a galaxy enters only if its stellar-kinematics file name matches
  an eCALIFA cube name (280 of the 300).
- **Common grid.** Stars and gas share the gas 0.5″ grid:
  - Its origin is the catalogue galaxy position, located through the cube
    WCS.
  - Each grid spaxel takes the values of the nearest 1″ stellar spaxel within
    0.75″.
  - The axes follow the MaNGA pipeline: x along increasing column, y along
    increasing row.

**Format check on excluded galaxies.** Readers are checked on galaxies that
are also in MaNGA. Those galaxies are excluded from C1, so inspecting them
costs no blinding. Their stellar and gas PAs must agree with MaNGA's to within
20°.
- CALIFA, UGC 08107 (MaNGA 11761-12705): stellar PA 300°, gas PA 309°,
  against MaNGA's 309° and 309°. This check is a unit test
  (`tests/test_c1.py`).
- SAMI, CATID 517164 (MaNGA 11754-1901): stellar PA 85°, gas PA 83°, against
  MaNGA's 83° and 83°. This check is also a unit test.

**Instrumental dispersion.** The pPXF stellar dispersions (both surveys) and
the LZIFU gas dispersions (SAMI; Zhou et al. 2017) are taken to be intrinsic,
i.e. instrument-corrected. Any residual error in this correction enters
counter-rotators and their matched controls alike. Because it appears only
through the k σ² term, it largely cancels in s.

**Overlap exclusion:**
- Any SAMI or CALIFA galaxy within 5″ of any MaNGA DR17 target is removed.
  The radius is 5″ rather than 3″ because CALIFA–MaNGA catalogue offsets
  reach 3″. This removes 101 galaxies from the parent sample.
- Galaxies in both SAMI and CALIFA (within 5″) are kept once, from CALIFA,
  because its wider field covers 1.5 R_e more often.

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

The sample is not chosen: it is every eligible galaxy in the public SAMI DR3
and CALIFA releases. No larger public integral-field sample independent of
MaNGA exists at present. At the expected N, a decisive result is about as
likely as not. A wrong
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
  q0 = 0.2; 30° ≤ i ≤ 80° required, as in the MaNGA run;
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

**Ring validity:** at least 20 valid spaxels of 0.5″ (5 arcsec²), covering
at least 3 azimuthal quadrants. This is the MaNGA rule; all maps are on
0.5″ grids.

**σ_*(1 R_e)** is the median valid stellar dispersion inside the deprojected
1 R_e ellipse, requiring at least 5 spaxels.

**Kinematic PAs** use every second spaxel in x and y, with at least 50 valid
spaxels, as in MaNGA.

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
- The evidence is a Bayes factor between two point hypotheses, H_V (a = P)
  and H_C (a = 0). With equal prior weight on the two points and a Gaussian
  likelihood for the intercept, it reduces to the likelihood ratio

  BF = exp[(a² − (a − P)²) / (2 σ_a²)],

  using the bootstrap σ_a as a plug-in standard error.

- The point estimate and 95% bootstrap percentile interval of f = a / P are also
  reported.

**Matching covariates** (standardised):
- log M_*;
- morphology (0–3 scale);
- log σ_*(1 R_e);
- inclination;
- survey (exact match: SAMI CR are matched to SAMI controls, and CALIFA CR to
  CALIFA controls).
- **Fallback rule.** If a survey's co-rotating pool is smaller than 3 × its CR
  count, every CR in that survey gets 2 controls. If the pool is smaller than
  2 × its CR count, every CR gets 1 control. Only if the pool is smaller than
  its CR count are that survey's CR excluded. The ratio used is reported.
- s_i uses the median D over those of its controls that have a valid ring
  at that radius. If none does, s_i is missing.

**Secondary analyses.** These are pre-specified and reported regardless of
outcome, but none of them decides the verdict:
1. As in the MaNGA run, S = median(D_CO) − median(D_CR) at 1.0 and 1.5 R_e,
   with bootstrap SD, at k−1, k and k+1.
2. b > 0, one-sided. This tests the conventional prediction that slowdown rises
   with twist. The test uses the one-sided 95% lower bound on b, i.e. the 5th
   percentile of the bootstrap distribution of b. The prediction holds if that
   bound is above 0.
3. The primary BF and f at k−1 and k+1.
4. The primary analysis for each survey separately.
5. A pooled MaNGA + SAMI + CALIFA analysis. It is labelled as combining
   exploratory and confirmatory data.
6. The Spearman correlation of s_i with predicted Rζ_i. H_V predicts a
   positive value.
7. The agreement of the pipeline's SAMI CR classification with the SAMI DR3
   catalogue PAs (PA_STELKIN, PA_GASKIN): the primary analysis is repeated
   on the CALIFA CR plus the SAMI CR whose published ΔPA also exceeds 150°.
8. The MaNGA exploratory primary statistic (intercept a, BF, f), recomputed
   without the Galaxy Zoo merger cut. This shows the size of the effect of
   omitting the merger cut.
9. Balance between CR and their matched controls. The table reports means and
   standardised mean differences for:
   - the matching covariates;
   - the stellar circular speed V_c,*(1 R_e);
   - the gas dispersion at 1.5 R_e.

   This checks that the drift correction acts similarly in both groups. It is
   reported, not used to modify the analysis.

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
  - are within 5″ of a MaNGA DR17 target;
  - lack either velocity map (SAMI files must match uniquely);
  - have a PA error above 20°;
  - have i < 30° or i > 80°;
  - have V_*(1 R_e) ≤ 40 km/s, or no valid stellar ring at 1 R_e;
  - lack a matching covariate (including SAMI TYPE 5 or −9, and CALIFA
    types outside the morphology map);
  - fail the catalogue-level quality flags listed under "Data collection
    procedures".
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
- Repository: https://github.com/JustinSchneider/vortex-chirality.
- Frozen commit: [[hash]], tagged `c1-prereg`.
- Archive: [[Zenodo DOI]].
- Entry points:
  - `analysis/c1_test.py select`: PAs only.
  - `analysis/c1_test.py dry`: counts only.
  - `analysis/c1_test.py run`: the registered analysis, run once.
  - `analysis/c1_secondary.py {5,7,8}`: secondary analyses.
- Supporting code: readers `src/sami.py` and `src/califa.py`; measurement
  `src/manga_measure.py` (unchanged from the MaNGA run).
- Any deviation from this plan will be listed, dated and justified in
  `preregistration/PREREG.md` and in the paper.

**Disclosure.** The analysis code and this plan were developed with an AI coding
assistant (Claude, Anthropic). All decisions are the author's.
