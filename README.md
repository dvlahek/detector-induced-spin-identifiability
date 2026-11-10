# Detector-restricted spin identifiability

Reproducible Python code and numerical datasets for studying the conditional identifiability of charged spin-1/2 and spin-1 particles from calibrated electromagnetic detector responses.

[![Reproduction CI](https://github.com/dvlahek/detector-induced-spin-identifiability/actions/workflows/reproduce.yml/badge.svg)](https://github.com/dvlahek/detector-induced-spin-identifiability/actions/workflows/reproduce.yml)
[![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

The primary model uses a one-dimensional signed-momentum detector response with effective on-shell measure $dk/(2E_k)$, without the $k^2$ Jacobian of a three-dimensional angular integral. It combines magnetic and longitudinal-current charge responses with no additional charge acceptance, $A_C(k)=1$. All bounds depend on the stated detector model and anchored Lipschitz form-factor class. This repository contains **source code and numerical evidence, not article or submission files**.

## Numerical results

| Check | Reference result |
| --- | --- |
| 99-point finite-bandwidth grid at $\delta=0$ | **63** sufficient whole-class separation certificates, **33** constructed overlaps, **3** unresolved |
| Moment-deformation scan with $|\delta|\le0.10$ | **61** sufficient whole-class certificates |
| Reference unit-covariance class distance | $0.4034833471\le d_{\rm class}\le0.4034946602$ |
| Global eight-segment numerical audit | 21 QCQPs, 12 positive numerical SOS lower bounds |
| High-resolution search | Four 48-segment profiles, including one known-overlap control and three unresolved cases |

The all-finite-bandwidth magnetic margin is strictly positive for a nonblind magnetic-moment domain and any finite Gaussian width and Lipschitz bound. See the [results and evidence map](docs/RESULTS.md) for assumptions and the exact scripts supporting each claim.

## Quick start

Python 3.12 is recommended. Run from the repository root on Linux, macOS, or Windows.

```bash
python -m pip install -r requirements.txt
python reproduce.py --mode smoke
python reproduce.py --mode full
```

- `smoke` checks frozen result counts and independently integrates the saved 48-segment profiles. The separate interval audit verifies the signs of all 33 constructive overlap endpoints.
- `full` regenerates the primary numerical tables, figure data, QCQP/SOS checks and estimator results.
- `highres` repeats the optional seeded 48-segment search with 200 starts per case.

```bash
python reproduce.py --mode highres
python theory_numerics/verify_overlap_sign_intervals.py
```

Generated results are written under `theory_numerics/reproduced_alpha0/`. The committed reference data remain unchanged. A positive result from local optimization is not proof of whole-class separation. The three unresolved cases remain unresolved.

## Repository contents

| Path | Purpose |
| --- | --- |
| [`reproduce.py`](reproduce.py) | Cross-platform entry point for verification and reproduction |
| [`theory_numerics/`](theory_numerics/) | Numerical implementations and validation scripts |
| [`theory_numerics/results/alpha0/`](theory_numerics/results/alpha0/) | Immutable reference tables, SOS witnesses, convergence data and 48-segment profiles |
| [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) | Installation, execution and verification |
| [`docs/RESULTS.md`](docs/RESULTS.md) | Claim-to-script-to-data evidence map |
| [`docs/SCOPE.md`](docs/SCOPE.md) | Assumptions and limitations |
| [`collider/`](collider/) and [`generator_configuration/`](generator_configuration/) | Auxiliary count-level benchmark and provenance |
| [`archive/development_snapshots/`](archive/development_snapshots/) | Historical numerical scripts and results, not used for the primary reproduction |

## Citation and scope

For fixed numerical results use the source revision described in [`CITATION.cff`](CITATION.cff), rather than a changing `main` branch. [Release and provenance information](docs/RELEASES.md) describes the frozen data and versions.

The results do not establish an operational collider detector, identify the spin of an observed particle, or settle the three unresolved grid cells. Numerical SOS feasibility is checked in floating point rather than verified interval arithmetic. See [scientific scope](docs/SCOPE.md) for details.

**Maintainer:** Dino Vlahek, Faculty of Organization and Informatics, University of Zagreb. Source code is released under the MIT license. No archival DOI is asserted.
