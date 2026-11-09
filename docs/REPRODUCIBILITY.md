# Reproducibility guide

## Frozen manuscript analysis

**Version:** `1.8.0`. **Exact scientific snapshot:** [`5c48258f69758818e3ed524cb4d7c1e9d5697811`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/5c48258f69758818e3ed524cb4d7c1e9d5697811). Use this SHA for the exact analysis cited in the manuscript, not the moving default branch.

```bash
git clone https://github.com/dvlahek/detector-induced-spin-identifiability.git
cd detector-induced-spin-identifiability
git checkout 5c48258f69758818e3ed524cb4d7c1e9d5697811
python3.12 -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

The CI workflow uses Python 3.12. Dependencies have **bounded version ranges, not a fully locked environment**; record `python --version` and `python -m pip freeze` for a portability audit.

## Verify the stored publication results

```bash
python theory_numerics/verify_publication_tables.py
```

This verifies committed numerical values; it is **not** a complete regeneration.

## Regenerate active theory outputs

Run from the repository root:

```bash
python theory_numerics/narrow_prior_magnetic_audit.py \
  --output reproduced_results/narrow_prior_magnetic_margins.csv

python theory_numerics/finite_bandwidth.py \
  --old-status theory_numerics/results/S1_slope_bandwidth_certificates_monotone_v1_1.csv \
  --outdir reproduced_results/publication_v1_0 \
  --unresolved-overlap-starts 1024

python theory_numerics/delta_robustness_audit.py \
  --status theory_numerics/results/publication_v1_0/finite_bandwidth_status.csv \
  --magnetic-reference theory_numerics/results/narrow_prior_magnetic_margins.csv \
  --outdir reproduced_results/delta_robustness

python theory_numerics/detector_measures_distance.py \
  --outdir reproduced_results/publication_v1_0 \
  --figdir reproduced_results/figures

python theory_numerics/verify_publication_tables.py \
  --data reproduced_results/publication_v1_0

python theory_numerics/precision_recovery.py \
  --output reproduced_results/S3_precision_based_recovery.csv
```

Static nuisance sensitivity can be recalculated via `python theory_numerics/nuisance_domain_sensitivity.py`. The main commands write to the ignored `reproduced_results/` folder, but this nuisance script updates its designated tracked results CSV.

## v1.8 response architecture and witnesses

```bash
python theory_numerics/charge_kernel_identity_check.py
python theory_numerics/explicit_overlap_family_verify.py
python theory_numerics/acceptance_sensitivity_audit.py
```

These scripts generate files beside the scripts under `theory_numerics/`. Use a disposable checkout if you need a pristine working tree. They verify the chosen positive acceptance, the 33 sign-changing two-kink witnesses and four **fixed** alternative acceptances. Details and evidence links: [RESULTS.md](RESULTS.md).

## CI and auxiliary material

The [reproduction workflow](../.github/workflows/reproduce.yml) runs the publication table checks, magnetic margins, 99-point classification, $\delta$, distance and v1.8 acceptance/witness audits, and estimator-space checks. It also runs **archival collider count-level checkpoints**; these are auxiliary and **not a validation input to the current theoretical manuscript**.

The public repository does not contain all collider event-level inputs or the modified UFO. Reproduction of the main theoretical paper requires neither of those. See [SCOPE.md](SCOPE.md) and [DATA_AVAILABILITY.md](../DATA_AVAILABILITY.md).

When reporting differences, give the exact commit, Python/dependency versions, invocation and expected/observed output. See [CONTRIBUTING.md](../CONTRIBUTING.md).
