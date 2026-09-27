# Detector-Restricted Identifiability of Elementary-Particle Spin

Code and processed data for the manuscript *Detector-Restricted Identifiability of Elementary-Particle Spin: Whole-Class Certificates and Exact Degeneracies*.

## Scientific question

An elementary particle has a definite spin, but an experiment observes responses of a specified interaction and detector. The calculation asks which spin distinctions remain identifiable when the electromagnetic form factors are uncertain and the detector has finite momentum bandwidth. The detailed example compares charged spin-1/2 and spin-1 candidates in a calibrated magnetic and charge-sensitive response model.

The absence of an internal rank-2 sector for spin 1/2 and the elementary-field constraint \(g_M+g_Q=1\) are theoretical inputs. The results here concern the detector resources and finite-bandwidth separation that follow from those inputs.

## Principal results

For the compact magnetic prior \(K=[1,3]\), the magnetic response has a strictly positive lower bound for every finite Gaussian bandwidth \(\sigma_k/m>0\) and finite Lipschitz constant \(L\). The accompanying 99-point calculation quantifies the margin over the bandwidth–smoothness grid; the smallest tabulated bound is approximately 0.005136 at \((0.50,24)\). The analytic statement is not limited to that grid.

For the broader \(g_M\in[-4,4]\) class, the same 99-point grid contains 64 whole-class separation certificates, 33 explicit origin overlaps and two unresolved points. The two unresolved cases are \((0.08,24)\) and \((0.10,16)\). The positive certificates use magnetic positivity and a pointwise charge envelope. Explicit overlaps use admissible continuous form factors. Finite-basis QCQP, SOS and multistart searches are recorded separately from these whole-class conclusions.

At \((\sigma_k/m,L)=(0.10,2)\), the unit-covariance whole-class distance is enclosed by \(0.4042553595\le d_{\rm class}\le0.4042657175\).

## Reproduction

Python 3.12 is recommended. From the repository root:

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

The default branch is the current code and data location cited by the manuscript. Authorship and software citation metadata are provided in [CITATION.cff](CITATION.cff). The code is released under the MIT license.
