# vortex-chirality

A chirality test of velocity-dependent explanations of galaxy rotation curves. Examples of such explanations are vorticity of space, gravitomagnetic-like fields and Coriolis-type terms.

A force that depends on the direction of motion treats stars and gas orbiting in opposite directions differently. The predicted asymmetry is ΔV = Rζ. Dark matter, MOND and modified inertia instead predict ΔV = 0 after the asymmetric-drift correction. Galaxies whose gas counter-rotates with respect to their stars test this directly.

## Status

- **Exploratory: MaNGA DR17.** 90 counter-rotators against 270 matched controls. The result is inconclusive as registered. Post hoc analysis shows that the counter-rotator "slowdown" follows the kinematic twist of the gas, not the predicted vortex amplitude.
- **Confirmatory: SAMI DR3 + CALIFA (C1).** Being pre-registered on OSF. See [preregistration/OSF_C1_confirmatory.md](preregistration/OSF_C1_confirmatory.md).

[preregistration/PREREG.md](preregistration/PREREG.md) is the append-only, dated analysis log. It records every amendment, the PA-convention bug that invalidated the first MaNGA run, and the corrected results. These dates are self-recorded and not independently timestamped. The C1 plan is timestamped by its OSF registration.

## Layout

- **derivations/**: sympy derivations, run as part of the tests.
  - `01_force_laws.py`: force laws and the chirality signature.
  - `02_dimensional_analysis.py`: Buckingham Π groups.
  - `03_magnitude_gates.py`: order-of-magnitude checks (seed, cosmic spin, coupling, lunar laser ranging, dark energy).
- **src/**: shared code.
  - `data.py`: SPARC loaders.
  - `models.py`: rotation-curve models and the chirality estimator.
  - `manga.py`: MaNGA sample and map access.
  - `manga_measure.py`: ring fits, drift correction and the D statistic.
- **analysis/**: one script per test. Each writes its output to `results/`.
  - `g1_manga_select.py` and `g1_manga_test.py`: the MaNGA test.
  - `g1_manga_robustness.py`, `g1_manga_discriminant_exploratory.py` and `s_pred_rar.py`: post hoc checks.
  - `c1_power.py`: power analysis for the confirmatory design.
- **data/raw/**: SPARC tables (Lelli, McGaugh & Schombert 2016) and the author's earlier SPARC rotation-curve fit tables.
- **data/processed/**: derived tables. MaNGA maps are not redistributed; they are fetched from the SDSS Science Archive Server.
- **docs/theory.md**: hypotheses, derivations and the verdict for each test.
- **manuscript/**: the paper, in preparation.

## Reproduce

```
python -m pytest                              # unit tests + derivations
python analysis/step0_foundation.py           # data/processed/sparc_canonical.csv
python analysis/g1_manga_select.py            # kinematic PAs (downloads MaNGA DAP maps)
python analysis/g1_manga_test.py              # the registered MaNGA test
python analysis/g1_manga_discriminant_exploratory.py
python analysis/c1_power.py
```

Dependencies: numpy, scipy, pandas, astropy, sympy, pafit, requests, pytest.

Analysis code and text were developed with an AI coding assistant (Claude, Anthropic).

## Licence

Code: MIT (see [LICENSE](LICENSE)). Third-party survey data remain under their own terms.
