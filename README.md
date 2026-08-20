# Detector-Induced Identifiability of Fundamental Spin

Public code, processed numerical data, sufficient statistics, and publication-facing provenance for the manuscript **Detector-Induced Identifiability of Fundamental Spin**.

This repository is intentionally separated from the private MadGraph/PYTHIA/Delphes development repository. The private repository was used to build and run the validated full-chain benchmark. Everything needed to verify the publication-facing numerical claims is collected here: finite-bandwidth validation code and tables, collider sufficient statistics, validated outputs, detector-readout analysis, exact process definitions, software versions, and reconstruction provenance.

## Version

Current publication snapshot: **v1.0.0**.

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
  figures/        finite-bandwidth status map
  monotone_closure.py
  verify_publication_tables.py
  plot_status_map.py

collider/
  sufficient_statistics_v2.json.gz.b64
  scripts/        publication analysis and public reproducer
  results/        validated numerical outputs
  figures/        publication-facing collider figures

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

python - <<'PY'
import base64, gzip
from pathlib import Path
src = Path('collider/sufficient_statistics_v2.json.gz.b64').read_text().strip()
Path('collider/sufficient_statistics_v2_decoded.json').write_bytes(
    gzip.decompress(base64.b64decode(src))
)
PY

python collider/scripts/reproduce_public_statistics.py \
  --statistics collider/sufficient_statistics_v2_decoded.json \
  --archived-closure collider/results/fermion_vector_kl_validation_summary.json \
  --archived-resource collider/results/architecture_kl_summary.csv \
  --output reproduced_results/public_statistics.json
```

The GitHub Actions workflow runs the same checks automatically. It verifies the finite-bandwidth `46/32/21` classification, decodes the lossless collider sufficient statistics, recomputes the public collider checkpoints, syntax-checks the complete event-level analysis scripts, and uploads the reproduced outputs.

The public collider reproducer verifies deterministic histogram divergences exactly and performs a fresh finite-sample bias/bootstrap calibration from the public counts. Exact frozen publication values are additionally checked against the archived validated outputs. The out-of-fold classifier cross-check needs event-level observables; its validated fold results are archived in `collider/results/fermion_vector_kl_validation_summary.json` rather than recomputed by the public CI.

## Finite-bandwidth archive

`theory_numerics/results/` contains the publication-facing S1-S7 tables. The final monotone 99-point status classification is explicit and machine-checkable. Because the admissible Lipschitz classes are nested in `L`, an explicit counterexample found at `L0` remains admissible for every `L >= L0` at fixed detector bandwidth. A failed counterexample search is never promoted to a proof of separation.

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

See `CITATION.cff`. The versioned GitHub Release `v1.0.0` is the publication snapshot; no external DOI archive is required for the manuscript's reproducibility statement.
