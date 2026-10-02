# VORTEX: a physical theory for the Rational Taper, and tests designed to break it

Start with [docs/theory.md](docs/theory.md), which covers the hypotheses, derivations, scaling law and results. All predictions and pass/fail thresholds were fixed in [preregistration/PREREG.md](preregistration/PREREG.md) before any data were analysed.

## Layout

- **derivations/** — sympy derivations:
  - `01_force_laws.py`: the force laws of both formulations
  - `02_dimensional_analysis.py`: the Buckingham Π analysis
  - `03_magnitude_gates.py`: order-of-magnitude checks
- **src/** — `data.py` (SPARC and RT-fit loaders) and `models.py` (the derived rotation-curve models)
- **analysis/** — one script per test; each writes its output to `results/*.json`
- **data/raw/** — SPARC tables and the RT fit tables, copied unchanged from `RT/rational-taper-validation`
- **data/counter_rotators.csv** — literature kinematics for Gate G1, with sources

## Reproduce

```
python -m pytest                        # model unit tests + all derivations
python analysis/step0_foundation.py     # builds data/processed/sparc_canonical.csv
python analysis/s2_scaling_law.py
python analysis/g4_flowing_space_fit.py
python analysis/g1_counter_rotation.py
```
