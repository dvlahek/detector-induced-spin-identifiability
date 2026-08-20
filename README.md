# Detector-Induced Identifiability of Fundamental Spin

Public code, numerical outputs, sufficient statistics, and publication-facing reproducibility material for the manuscript **Detector-Induced Identifiability of Fundamental Spin**.

This repository is intentionally separated from the private MadGraph/PYTHIA/Delphes development repository. The private repository was used to produce the validated full-chain benchmark. This public repository contains the final analysis code, exact generator/detector provenance, publication outputs, and histogram sufficient statistics needed to verify the primary collider and detector-readout claims without exposing the development repository.

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
- recentered-bootstrap 95% interval: `[0.1447812963, 0.1529124751]`;
- classifier cross-check: `0.1736703 +/- 0.0040008`.

Detector-readout ablation at the same geometric acceptance:

- angles + energy: corrected conditional minimum KL about `0.1613`;
- two directions: about `0.00468`;
- one direction: about `0.00116`;
- energy only: about `0.1527`.

The collider readout ablation is a separate detector-architecture study. It is **not** identified with the primitive magnetic and charge-sensitive resources in the analytic protected-sector theorem.

## Repository layout

```text
theory_numerics/
  results/        finite-bandwidth publication tables
  figures/        reproduced finite-bandwidth status map
  monotone_closure.py
  verify_publication_tables.py
  plot_status_map.py

collider/
  sufficient_statistics.json
  full3d_12_counts.json
  scripts/        publication analysis and public reproducer
  results/        validated numerical outputs

generator_configuration/
  final process definitions, PYTHIA command files, model metadata,
  software versions, detector-card provenance, and workflow provenance
```

## Reproduce the public checkpoints

Python 3.12 is recommended.

```bash
python -m pip install -r requirements.txt
python theory_numerics/verify_publication_tables.py
python theory_numerics/plot_status_map.py
python collider/scripts/reproduce_public_statistics.py \
  --statistics collider/sufficient_statistics.json \
  --full12 collider/full3d_12_counts.json \
  --archived-closure collider/results/fermion_vector_kl_validation_summary.json \
  --archived-resource collider/results/architecture_kl_summary.csv \
  --output reproduced_results/public_statistics.json
```

The GitHub Actions workflow runs the same checks automatically.

The public collider reproducer verifies the deterministic raw histogram divergences exactly and performs a fresh finite-sample bias/bootstrap calibration from the public counts. The exact frozen publication values are additionally checked against the archived validated outputs. The out-of-fold classifier cross-check needs event-level observables, so its validated fold results are archived in `collider/results/fermion_vector_kl_validation_summary.json` rather than rerun by the public CI.

## Finite-bandwidth archive

`theory_numerics/results/` contains the publication-facing S1--S7 tables. The final monotone 99-point status classification is explicit and machine-checkable. Because the admissible Lipschitz classes are nested in `L`, an explicit counterexample found at `L0` remains admissible for every `L >= L0` at fixed detector bandwidth. A failed counterexample search is never promoted to a proof of separation.

## Generator and detector provenance

The validated private production workflow used:

- MadGraph5_aMC@NLO 3.7.2 pinned commit `be1e7b273ca961c335ff2ee6da3688b5049b069e`;
- PYTHIA 8.312;
- ROOT 6.40.02;
- Delphes 3.5.1 with the CLICdet Stage-1 card;
- 100k closure run ID `31995375143`;
- source commit `58210d4e7596a2752a3512dcf7908372d287b736`;
- detector-readout workflow run ID `32001628799`;
- readout-analysis commit `cbe53a8cfeb3ec22086a43baf8430b2297e739cf`.

The exact final process definitions and PYTHIA command files are in `generator_configuration/`. The private development repository itself is not required to verify the publication-facing numerical claims.

## Citation

See `CITATION.cff`. A versioned archival DOI can be added to the manuscript Data Availability Statement after the submission release is archived.
