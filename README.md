# Whole-Class Certificates for Detector-Restricted Identifiability of Elementary-Particle Spin

Clean publication 1.0 reproducibility repository for the PRD-focused manuscript.
It contains the active whole-class calculations, machine-readable certificates,
collider sufficient statistics, and generator/reconstruction provenance.

## Publication baseline

The canonical publication baseline is **1.0.0** (`VERSION = 1.0.0`) on the
default `main` branch:

https://github.com/dvlahek/detector-induced-spin-identifiability

Superseded development material is retained under
`archive/development_snapshots/`; the repository root exposes one current
publication path.

## Main checkpoints

- 99 finite-bandwidth points: 64 whole-class certificates, 33 explicit
  overlaps, and 2 unresolved points.
- Magnetic positivity forces every explicit protected-origin overlap to the
  broad-class blind value `g_M=0`; none is an overlap at the natural `g_M=2`
  vector point.
- At `(sigma_k/m,L)=(0.10,2)`, the unit-covariance class distance is tightly
  enclosed by `0.4042553595 <= d_class <= 0.4042657175`.
- The rigorous lower distance bound is tabulated on all 99 points and is
  positive at exactly the same 64 certified points.
- The explicit 32-segment overlap at `(0.15,8)` is admissible in both the
  covariant and Sachs slope conventions: `max|G_C'|=7.927785 < 8`.
- The same 1024-start, 32-segment search remains positive at `(0.08,24)` and
  `(0.10,16)`; these remain search outcomes, not proofs.
- The matched detector-level reference gives bias-corrected minimum KL
  `0.1488976547`, with bootstrap interval `[0.1447812963,0.1529124751]`.

The collider calculation is a fixed-energy detector-readout closure, not a
gauge-complete vector theory and not a test of the analytic theorem.

## Active layout

```text
theory_numerics/
  finite_bandwidth.py
  detector_measures_distance.py
  precision_recovery.py
  verify_publication_tables.py
  nuisance_domain_sensitivity.py
  results/publication_v1_0/

collider/
  sufficient_statistics_v2.json.gz.b64
  scripts/
  results/

generator_configuration/
  exact process cards, software versions, seeds, and provenance

archive/development_snapshots/
  legacy workflows, release notes, scripts, and versioned numerical outputs
```

## Reproduce publication 1.0

Python 3.12 is recommended.

```bash
python -m pip install -r requirements.txt

python theory_numerics/nuisance_domain_sensitivity.py

python theory_numerics/finite_bandwidth.py \
  --old-status theory_numerics/results/S1_slope_bandwidth_certificates_monotone_v1_1.csv \
  --outdir reproduced_results/publication_v1_0 \
  --unresolved-overlap-starts 1024

python theory_numerics/detector_measures_distance.py \
  --outdir reproduced_results/publication_v1_0 \
  --figdir reproduced_results/figures

python theory_numerics/verify_publication_tables.py \
  --data reproduced_results/publication_v1_0

python theory_numerics/precision_recovery.py \
  --output reproduced_results/S3_precision_based_recovery.csv
```

The GitHub Actions workflow also decodes the lossless public collider
sufficient statistics, recomputes the publication-facing collider checkpoints,
and syntax-checks the event-level analyses.

## Citation

The repository URL above is the canonical code and data location. See
`CITATION.cff` for authorship metadata.
