# Detector-restricted spin identifiability

Research code, numerical data and LaTeX sources for *Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles: Whole-Class Certificates and Exact Degeneracies*.

[![Reproduction CI](https://github.com/dvlahek/detector-induced-spin-identifiability/actions/workflows/reproduce.yml/badge.svg)](https://github.com/dvlahek/detector-induced-spin-identifiability/actions/workflows/reproduce.yml) [![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Status:** Manuscript prepared for *Journal of Physics A*. Not yet peer reviewed or published.

**Reference model:** the positive magnetic and longitudinal-current detector responses, with no additional momentum-dependent charge acceptance ($A_C=1$). All numerical claims are conditional on the detector model and specified anchored Lipschitz form-factor class.

## Main results

| Result | Finding |
| --- | --- |
| Nonblind magnetic-moment class | Separation at every finite Gaussian bandwidth and finite Lipschitz bound |
| Broad 99-point grid ($\delta=0$) | **63** sufficient whole-class certificates, **33** constructed overlaps, **3** unresolved |
| Moment deformation $|\delta|\le0.10$ | **61** sufficient whole-class certificates |
| Reference unit-covariance distance | $0.4034833471\le d_{\rm class}\le0.4034946602$ |
| Independent numerical checks | 21 eight-segment QCQPs, 12 positive numerical SOS bounds and a 48-segment adversarial search |

![99-point classification](manuscript_jphysa/figures/figure2a.png)

Green indicates a sufficient whole-class separation certificate, red a constructed overlap and grey a case left unresolved by the tests.

## Reproduce the results

Python 3.12 is recommended. From the repository root on Windows, Linux or macOS:

```bash
python -m pip install -r requirements.txt
python reproduce.py --mode smoke
python reproduce.py --mode full
```

The **smoke** check validates frozen result counts and independently re-integrates four saved 48-segment form-factor profiles. The **full** run regenerates the primary tables, figures, estimator results and numerical certificate checks. Generated outputs go to `theory_numerics/reproduced_alpha0/` and do not overwrite the published data.

The adversarial search with 200 seeded starts per case can also be run independently:

```bash
python reproduce.py --mode highres
```

A positive optimized response does **not** prove separation. The three unresolved cases remain unresolved.

Full setup, expected outputs and a fresh LaTeX build are documented in [Reproducibility](docs/REPRODUCIBILITY.md).

## Find the material

| Location | Contents |
| --- | --- |
| [`manuscript_jphysa/`](manuscript_jphysa/) | Main text, Supplementary Material, bibliography and figures |
| [`theory_numerics/spin_alpha0_*.py`](theory_numerics/) | Primary finite-width computations |
| [`theory_numerics/results/alpha0/`](theory_numerics/results/alpha0/) | Frozen tables, SOS witnesses, 48-segment profiles and convergence data |
| [`reproduce.py`](reproduce.py) | Cross-platform command-line entry point |
| [Evidence map](docs/RESULTS.md) | Claim-to-code-to-data references |
| [Scientific scope](docs/SCOPE.md) | Assumptions, proof status and limitations |

## Cite the software

The fixed, reproducible scientific snapshot is [`93c70ed48244`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/93c70ed482443730c1abc6c2891732fb1bd5bd4f). Use [`CITATION.cff`](CITATION.cff) to cite the code and data. Do not substitute the changing `main` branch for the frozen snapshot.

**Author:** Dino Vlahek, Faculty of Organization and Informatics, University of Zagreb.

The reported results do not demonstrate an actual collider detector, determine the spin of an observed particle or settle the three unresolved cells. The license applies to the code. No Zenodo DOI is claimed.
