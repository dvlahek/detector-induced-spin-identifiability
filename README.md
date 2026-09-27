# Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles

Code and processed data for the manuscript *Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles: Whole-Class Certificates and Exact Degeneracies*.

## Problem and model

An experiment measures currents, polarization and reconstructed particles, not a spin label. Here we determine which electromagnetic responses distinguish charged spin-1/2 and spin-1 candidates after detector bandwidth and uncertain form factors are taken into account. The calculation concerns the calibrated magnetic and charge-sensitive response model defined in the manuscript.

The missing internal rank-2 sector of spin 1/2 and Haberzettl's elementary-field relation \(g_M+g_Q=1\) are theoretical inputs. The contribution here is the detector-resource criterion and its finite-bandwidth analysis over the admitted form-factor class.

## Results

For the compact magnetic prior \(K=[1,3]\), the magnetic response has a strictly positive lower bound for every finite Gaussian bandwidth \(\sigma_k/m>0\) and finite Lipschitz constant \(L\). The accompanying 99-point calculation quantifies the margin over the bandwidth–smoothness grid; the smallest tabulated bound is approximately 0.005136 at \((0.50,24)\). The analytic statement is not limited to that grid.

For the broader \(g_M\in[-4,4]\) class, the same 99-point grid contains 64 whole-class separation certificates, 33 explicit origin overlaps and two unresolved points. The two unresolved cases are \((0.08,24)\) and \((0.10,16)\). The positive certificates use magnetic positivity and a pointwise charge envelope. Explicit overlaps use admissible continuous form factors. Finite-basis QCQP, SOS and multistart searches are recorded separately from these whole-class conclusions.

At \((\sigma_k/m,L)=(0.10,2)\), the unit-covariance whole-class distance is enclosed by \(0.4042553595\le d_{\rm class}\le0.4042657175\).

## Reproduce the calculations

The manuscript cites code version `1.6.0` at commit [`80d3bfd8`](https://github.com/dvlahek/detector-induced-spin-identifiability/commit/80d3bfd8d87a4c011e651fdbd90345f22842441c). Use that commit for the exact submitted calculation. The older GitHub releases document earlier stages and are not the manuscript's computational version. Python 3.12 is recommended.

From the repository root:

```bash
python -m pip install -r requirements.txt

# Compact-prior magnetic margins
python theory_numerics/narrow_prior_magnetic_audit.py \
  --output reproduced_results/narrow_prior_magnetic_margins.csv

# Broad-class finite-bandwidth classification
python theory_numerics/finite_bandwidth.py \
  --old-status theory_numerics/results/S1_slope_bandwidth_certificates_monotone_v1_1.csv \
  --outdir reproduced_results/publication_v1_0 \
  --unresolved-overlap-starts 1024

# Detector measures, distance bounds and coordinate comparison
python theory_numerics/detector_measures_distance.py \
  --outdir reproduced_results/publication_v1_0 \
  --figdir reproduced_results/figures

python theory_numerics/verify_publication_tables.py \
  --data reproduced_results/publication_v1_0

# Static nuisance-domain sensitivity and estimator-space audit
python theory_numerics/nuisance_domain_sensitivity.py
python theory_numerics/precision_recovery.py \
  --output reproduced_results/S3_precision_based_recovery.csv
```

The GitHub Actions workflow checks the published numerical checkpoints, the 99 compact-prior margins and the public collider sufficient statistics. The generated files are written to `reproduced_results/`, except for the static nuisance-domain script, which updates its designated results CSV.

## Repository contents

- `theory_numerics/`: form-factor calculations, certified bounds, constructive overlaps, numerical audits and their processed results.
- `collider/`: fixed-energy detector-readout example, public sufficient statistics and reconstructed-level validation.
- `generator_configuration/`: process definitions, model metadata, software versions, detector card provenance and seeds.
- `archive/development_snapshots/`: earlier computational snapshots retained for provenance; current reproduction uses the active scripts above.

The collider example uses photon-channel production at one fixed energy. Its public count-level data reproduce the reported histogram diagnostics and readout ablations. The event-level sample and modified vector UFO are not distributed; the event-level classifier result is provided as a validated numerical record, not as a publicly rerunnable computation. This example does not implement the two primitive protected electromagnetic response blocks used in the analytic certificate.

## Citation and availability

The manuscript's fixed reference is the commit linked above. The default branch keeps the current code and data, including documentation updates. Software citation metadata are in [CITATION.cff](CITATION.cff), and the code is released under the MIT license.
