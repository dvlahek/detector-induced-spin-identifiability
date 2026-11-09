# Detector-restricted spin identifiability

**Research code and reproducible numerical results for the conditional identification of charged spin-1/2 and spin-1 particles from electromagnetic detector responses.**

[![Reproduction CI](https://github.com/dvlahek/detector-induced-spin-identifiability/actions/workflows/reproduce.yml/badge.svg)](https://github.com/dvlahek/detector-induced-spin-identifiability/actions/workflows/reproduce.yml)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Research snapshot v1.9.0](https://img.shields.io/badge/research%20snapshot-v1.9.0-informational.svg)](CITATION.cff)

**Manuscript:** *Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles: Whole-Class Certificates and Exact Degeneracies* — prepared for journal submission; **not yet published or peer-reviewed**.  
**Author:** Dino Vlahek (Faculty of Organization and Informatics, University of Zagreb).  
**Frozen scientific code/data:** [`46271df7ba4b`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/46271df7ba4b94cfa6ad2af7f2627209c1e0ae7b).  
**Primary charge channel:** unfiltered squared-current response, $A_C(k)=1$ ($\alpha=0$).

[Manuscript sources](manuscript_jphysa/) · [Evidence map](docs/RESULTS.md) · [Reproduce](docs/REPRODUCIBILITY.md) · [Scope](docs/SCOPE.md) · [Citation](#citation)

## Research question

Spin is intrinsic, but a detector observes only its available response blocks. We ask **which calibrated electromagnetic channels can distinguish charged spin-1/2 from spin-1 models when form factors are uncertain and detector bandwidth is finite**. The lower-spin internal rank-2 sector is absent; the higher-spin current is restricted by an explicit moment/nuisance class.

The finite-width *primary* charge measure is the normalized squared-longitudinal-current weight, $d\nu_C\propto\rho_0(r)r^2x\,dr$. There is no additional, momentum-dependent charge-event acceptance in the main model. The formerly primary $\alpha=2$ acceptance is retained **only as a sensitivity case**.

> **Scope.** The results are conditional on the stated positive detector-response model and anchored Lipschitz form-factor class. This repository does not establish a complete process-dependent detector implementation, a particle discovery, or a universal statement across all interactions. The collider folders record auxiliary historical experiments, not the numerical basis of this manuscript.

## Main results

- **All-finite-bandwidth magnetic separation:** if $\inf_{g\in K}|g|>0$, a strictly positive magnetic-only margin exists at every finite Gaussian bandwidth and finite Lipschitz constant, independently of $\delta$.
- **Broad-class 99-point grid** (primary $\alpha=0$): **63 whole-class sufficient certificates, 33 explicit continuous-class origin overlaps, 3 unresolved cases**. Unresolved does not imply degeneracy or separation.
- **Moment deformation:** for $|\delta|\leq0.10$, **61** of the same grid points have a uniform whole-class sufficient certificate.
- **Reference finite-width distance** at $(\sigma_k/m,L)=(0.10,2)$: $0.4034833471\leq d_{\mathrm{class}}\leq0.4034946602$ (unit covariance).
- **Fixed acceptance sensitivity** $A_C^{(\alpha)}=[2m/(E_k+m)]^\alpha$: $\alpha=0,1,2,3$ gives (63,33,3), (63,33,3), (64,33,2), (64,33,2), respectively. This is not a theorem uniform in $\alpha$.

![Unfiltered charge-response 99-point classification](manuscript_jphysa/figures/figure2a.png)

*Green: sufficient whole-class certificate; red: signed constructive overlap; gray: unresolved. See [assumptions and detailed evidence](docs/RESULTS.md).*

## Quick start

With Python 3.12, from repository root:

```bash
python -m pip install -r requirements.txt
python theory_numerics/spin_alpha0_compute.py
python theory_numerics/spin_alpha0_qcqp.py
python theory_numerics/spin_alpha0_power.py
python theory_numerics/verify_alpha0_release.py
```

The results are generated in the ignored `theory_numerics/reproduced_alpha0/` directory and compared against the fixed reference data under [`theory_numerics/results/alpha0/`](theory_numerics/results/alpha0/). See the [full reproducibility guide](docs/REPRODUCIBILITY.md).

## Repository map

| Path | Role in the current paper |
| --- | --- |
| [`manuscript_jphysa/`](manuscript_jphysa/) | Main/supplement LaTeX and figure assets |
| [`theory_numerics/spin_alpha0_*.py`](theory_numerics/) | Primary unfiltered calculations |
| [`theory_numerics/results/alpha0/`](theory_numerics/results/alpha0/) | Fixed numerical evidence and SOS witness records |
| [`docs/`](docs/) | Scientific scope, evidence mapping, reproduction, version history |
| [`theory_numerics/results/publication_v1_0/`](theory_numerics/results/publication_v1_0/) | **Legacy filtered-charge** snapshot, not current results |
| [`collider/`](collider/) and [`generator_configuration/`](generator_configuration/) | Auxiliary development provenance (not used in paper) |
| [`archive/development_snapshots/`](archive/development_snapshots/) | Historical development archive |

## Citation

Use [`CITATION.cff`](CITATION.cff) or GitHub's *Cite this repository* control, citing the immutable science commit rather than changing `main`:

> Vlahek, D. (2026). *Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles: Whole-Class Certificates and Exact Degeneracies — Code and Data*, v1.9.0, commit [`46271df7ba4b`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/46271df7ba4b94cfa6ad2af7f2627209c1e0ae7b), GitHub.

No DOI is asserted. The MIT license covers repository code; see [data availability](DATA_AVAILABILITY.md) and [contribution guidance](CONTRIBUTING.md).
