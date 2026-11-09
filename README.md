# Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles

Code and processed data for the manuscript *Detector-Restricted Identifiability of Charged Spin-1/2 and Spin-1 Particles: Whole-Class Certificates and Exact Degeneracies*.

## Problem and model

An experiment measures currents, polarization and reconstructed particles, not a spin label. Here we determine which electromagnetic responses distinguish charged spin-1/2 and spin-1 candidates after detector bandwidth and uncertain form factors are taken into account. The calculation concerns the calibrated magnetic and charge-sensitive response model defined in the manuscript.

The missing internal rank-2 sector of spin 1/2 and Haberzettl's elementary-field relation \(g_M+g_Q=1\) define the reference current class. The current analysis also tests the deformation \(g_M+g_Q=1+\delta\). The contribution is the detector-resource criterion and its finite-bandwidth analysis over complete admitted form-factor classes.

## Results

For the compact magnetic prior \(K=[1,3]\), the magnetic response has a strictly positive lower bound for every finite Gaussian bandwidth \(\sigma_k/m>0\) and finite Lipschitz constant \(L\). The accompanying 99-point calculation quantifies the margin over the bandwidth–smoothness grid; the smallest tabulated bound is approximately 0.005136 at \((0.50,24)\). The analytic statement is not limited to that grid.

For the broader \(g_M\in[-4,4]\) class, the same 99-point grid contains 64 whole-class separation certificates, 33 explicit origin overlaps and two unresolved points. The two unresolved cases are \((0.08,24)\) and \((0.10,16)\). The positive certificates use magnetic positivity and a pointwise charge envelope. Explicit overlaps use admissible continuous form factors. Finite-basis QCQP, SOS and multistart searches are recorded separately from these whole-class conclusions.

At \((\sigma_k/m,L)=(0.10,2)\), the \(\delta=0\) unit-covariance whole-class distance is enclosed by \(0.4042553595\le d_{\rm class}\le0.4042657175\). For the enlarged class \(|\delta|\le0.10\), 61 of the 99 points retain a whole-class certificate and the reference charge-envelope lower bound is 0.737365. The compact-prior magnetic theorem is independent of \(\delta\).

## Calibrated charge acceptance and robustness (v1.8.0)

The model uses a positive, channel-dependent acceptance A_C(k)=[2m/(E_k+m)]^2. This is a chosen detector architecture, not a universal consequence of the spin-one current. The two-kink witness family provides constructive origin overlaps in all 33 designated cells. A fixed-alpha sensitivity audit tests A_C^(alpha) for alpha=0,1,2,3 with certified/overlap/unresolved counts 63/33/3, 63/33/3, 64/33/2, 64/33/2. These are separate model checks, not a certificate uniform in alpha. The main mathematical result for the magnetic channel is acceptance-independent.

Reproduction files: `theory_numerics/explicit_overlap_family_verify.py`, `theory_numerics/acceptance_sensitivity_audit.py`, `theory_numerics/charge_kernel_identity_check.py` and processed CSV/JSON files under `theory_numerics/results/`. Full source of the manuscript and figures is under `manuscript_jphysa/` in the submission branch, with bibliography pointing to the immutable numerical-code commit.

## Reproduce the calculations

The current J. Phys. A computational version is `1.8.0`. The manuscript cites a fixed commit SHA so the submitted calculation remains reproducible even if the default branch later changes. Older GitHub releases document earlier stages and are not the current manuscript version. Python 3.12 is recommended.

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

# Uniform robustness to g_M + g_Q = 1 + delta
python theory_numerics/delta_robustness_audit.py \
  --status theory_numerics/results/publication_v1_0/finite_bandwidth_status.csv \
  --magnetic-reference theory_numerics/results/narrow_prior_magnetic_margins.csv \
  --outdir reproduced_results/delta_robustness

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

The GitHub Actions workflow checks the manuscript numerical checkpoints, the 99 compact-prior margins, and the \(\delta\)-robustness summary. The collider directory is retained as an auxiliary development record and is not part of the current J. Phys. A manuscript. The generated files are written to `reproduced_results/`, except for the static nuisance-domain script, which updates its designated results CSV.

## Repository contents

- `theory_numerics/`: form-factor calculations, certified bounds, constructive overlaps, numerical audits and their processed results.
- `collider/`: auxiliary fixed-energy readout study retained for provenance; it is not used by the current manuscript.
- `generator_configuration/`: provenance for that auxiliary collider study.
- `archive/development_snapshots/`: earlier computational snapshots retained for provenance; current reproduction uses the active scripts above.

## Data availability

The data used by the current manuscript are described in [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). The theoretical results are publicly reproducible from the processed data and scripts in `theory_numerics/`.

## Citation and availability

The manuscript's software reference records the fixed submission commit. The default branch keeps the current code and data, including documentation updates. Software citation metadata are in [CITATION.cff](CITATION.cff), and the code is released under the MIT license.


Additional J. Phys. A v1.8.0 checks:

```bash
python theory_numerics/charge_kernel_identity_check.py
python theory_numerics/explicit_overlap_family_verify.py
python theory_numerics/acceptance_sensitivity_audit.py
```
