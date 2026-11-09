# Reproduce the J. Phys. A numerical results (v1.9.0)

**Code/data SHA:** [`f711505b03b1a01ad4d4be8a50ce4d46dbda2080`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/f711505b03b1a01ad4d4be8a50ce4d46dbda2080). This revision uses the **unfiltered $\alpha=0$** charge measure. Older scripts under `theory_numerics/` and datasets under `results/publication_v1_0/` are filtered $\alpha=2$ provenance and are not part of the current numerical reproduction.

## Environment

Python 3.12 is recommended. Dependencies are range-pinned rather than fully lockfile-pinned; report exact versions for cross-platform comparisons.

```bash
git clone https://github.com/dvlahek/detector-induced-spin-identifiability.git
cd detector-induced-spin-identifiability
git checkout f711505b03b1a01ad4d4be8a50ce4d46dbda2080
python -m pip install -r requirements.txt
```

## Reproduce and compare the frozen primary results

```bash
python theory_numerics/spin_alpha0_compute.py
python theory_numerics/spin_alpha0_qcqp.py
python theory_numerics/spin_alpha0_power.py
python theory_numerics/spin_alpha0_figures.py
python theory_numerics/verify_alpha0_release.py
```

Outputs are written under ignored `theory_numerics/reproduced_alpha0/`. The verifier compares regenerated results with [`theory_numerics/results/alpha0/`](../theory_numerics/results/alpha0/), asserting 63/33/3, 46 componentwise certificates, 62 Sachs certificates, the seven $\Delta$ counts, 21 positive finite-basis QCQP minima and at least 12 positive numerical SOS feasibility checks.

The 33 witness signs are evaluated with knot-split quadrature; their exact roots follow by continuity if those signs hold. The SOS matrix data are numerical, not interval-arithmetic certificates. Figure generation requires Matplotlib; estimator power uses a fixed random seed, specified finite sampling and nearest-neighbor cloud evaluation.

## Auxiliary and archived calculations

The magnetic-only theorem is independently checked with:

```bash
python theory_numerics/narrow_prior_magnetic_audit.py \
  --output reproduced_results/narrow_prior_magnetic_margins.csv
```

The four fixed $\alpha$ detector-acceptance cases are available in [`theory_numerics/results/alpha0/acceptance_sensitivity_99point.csv`](../theory_numerics/results/alpha0/acceptance_sensitivity_99point.csv). The earlier $\alpha=2$ analysis remains available at commit `5c48258f69758818e3ed524cb4d7c1e9d5697811` for historical comparison, **not** as current evidence.
