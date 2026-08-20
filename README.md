# Detector-Induced Identifiability of Fundamental Spin

Public code, frozen numerical data, and publication-facing reproducibility material for the manuscript **Detector-Induced Identifiability of Fundamental Spin**.

This repository is intentionally separated from the private MadGraph/PYTHIA/Delphes development repository. The public analysis starts from a frozen reconstructed 100k-event reference sample and includes the exact validation and detector-readout analysis used for the manuscript, together with generator provenance sufficient to identify the production configuration.

## Main publication checkpoints

Finite-bandwidth protected-sector scan:

- 99 `(sigma_k/m, L)` points;
- 46 rigorously certified points;
- 32 explicit-overlap points after monotone closure of nested slope classes;
- 21 unresolved points.

Matched detector-level fermion-versus-vector closure:

- fermion acceptance: `0.92758`;
- vector acceptance: `0.91766`;
- 12^3 + lost raw minimum KL: `0.1559913222`;
- estimated finite-N plug-in bias: `0.0070936675`;
- bias-corrected minimum KL: `0.1488976547`;
- recentered bootstrap 95% interval: `[0.1447812963, 0.1529124751]`;
- classifier cross-check: `0.1736703 +/- 0.0040008`.

Frozen detector-readout ablation at the same geometric acceptance:

- angles + energy: corrected conditional minimum KL about `0.1613`;
- two directions: about `0.00468`;
- one direction: about `0.00116`;
- energy only: about `0.1527`.

The collider readout ablation is a separate detector-architecture study. It is **not** identified with the primitive magnetic and charge-sensitive resources in the analytic protected-sector theorem.

## Repository layout

```text
theory_numerics/
  results/        archived finite-bandwidth publication tables
  figures/        finite-bandwidth status map
  monotone_closure.py
  verify_publication_tables.py
  plot_status_map.py

collider/
  frozen_data/    frozen 100k reconstructed reference
  scripts/        exact publication analysis scripts
  results/        validated numerical outputs
  figures/        validation and readout-ablation figures

generator_configuration/
  exact process definitions, PYTHIA command files, model metadata,
  software versions, detector-card provenance, and workflow provenance
```

## Reproduce the publication-facing analysis

Python 3.12 is recommended.

```bash
python -m pip install -r requirements.txt
```

### 1. Verify finite-bandwidth publication tables

```bash
python theory_numerics/verify_publication_tables.py
python theory_numerics/plot_status_map.py
```

The archived monotone classification contains exactly `46 / 32 / 21` guaranteed / explicit-overlap / unresolved points.

### 2. Recompute the collider KL validation

```bash
python collider/scripts/validate_fermion_vector_kl.py \
  --reference collider/frozen_data/reco_reference_100k.npz \
  --outdir reproduced_results/closure \
  --binnings 6 8 10 12 \
  --replicates 300 \
  --bootstrap-replicates 500 \
  --classifier-max-per-class 100000 \
  --seed 20260825
```

### 3. Recompute the detector-readout degradation study

```bash
python collider/scripts/resource_degradation_analysis.py \
  --reference collider/frozen_data/reco_reference_100k.npz \
  --outdir reproduced_results/resource \
  --seed 20260817 \
  --bins 12 \
  --pseudo 0.5 \
  --bias-reps 500 \
  --bootstrap-reps 1000 \
  --mc-trials 50000
```

The resource-study Monte Carlo was generated on GitHub Actions with Python 3.12.13, NumPy 2.5.2, pandas 3.0.5, and Matplotlib 3.11.1. These versions are pinned because multinomial random streams can vary across NumPy versions even when the deterministic raw KL is unchanged.

## Finite-bandwidth numerical archive

`theory_numerics/results/` contains the final publication-facing S1--S7 tables, including the monotone 99-point status classification and the covariance, estimator-recovery, ablation, negative-control, basis-convergence, and Standard Model W-point summaries.

The status-map table and figure are the archived final monotone classification used by the manuscript. The monotone closure is explicit in `monotone_closure.py`: because the admissible Lipschitz classes are nested in `L`, an explicit counterexample found at `L0` remains admissible for every `L >= L0` at fixed detector bandwidth. A failed counterexample search is never promoted to a proof of separation.

## Frozen collider reference

`collider/frozen_data/reco_reference_100k.npz` contains the publication reference used by both collider analyses. It stores the three reconstructed observables for accepted events and the generated event counts. No total-rate information is needed for the readout-degradation study, which conditions on a fixed accepted event budget.

The event-generation infrastructure itself remains private. Its exact final configuration is documented in `generator_configuration/`; the public frozen sample allows the manuscript's statistical and readout analyses to be rerun without MadGraph, PYTHIA, ROOT, or Delphes.

## Provenance

The validated private production workflow used:

- 100k closure run ID `31995375143`;
- source commit `58210d4e7596a2752a3512dcf7908372d287b736`;
- detector-readout workflow run ID `32001628799`;
- readout-analysis commit `cbe53a8cfeb3ec22086a43baf8430b2297e739cf`.

Checksums for the frozen input and archived outputs are in `CHECKSUMS.sha256`.

## Citation

See `CITATION.cff`. A versioned archival DOI can be added to the manuscript Data Availability Statement after the submission release is archived.
