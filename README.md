# Whole-Class Certificates for Detector-Restricted Identifiability of Elementary-Particle Spin

Public code, processed numerical data, sufficient statistics, and provenance
for the PRD-focused manuscript **Whole-Class Certificates for
Detector-Restricted Identifiability of Elementary-Particle Spin**.

The public reproduction path is deliberately separated from the private
MadGraph/PYTHIA/Delphes development repository. Everything needed to verify the
publication-facing numerical claims is collected here: whole-class and
finite-basis certificates, finite-bandwidth tables, collider sufficient
statistics, validated detector-readout outputs, exact process definitions,
software versions, and reconstruction provenance.

## Version

Current fixed publication snapshot: **v1.5.0**.

Permanent release URL:
https://github.com/dvlahek/detector-induced-spin-identifiability/releases/tag/v1.5.0

## Main publication checkpoints

Static nuisance-domain sensitivity for the protected spin-1 fingerprint
`(W_M,W_C)=(g^2,1-g)`:

- magnetic-only uniform separation fails iff `0 in K`;
- charge-only uniform separation fails iff `1 in K`;
- the combined protected response remains separated for every nonempty compact
  `K`;
- the minimum primitive support is `tau(K)=1` if at least one blind point is
  excluded and `tau(K)=2` only when both `0` and `1` are admitted;
- for `K_R=[2-R,2+R]`, the transition to two required primitive resources is
  exactly at `R=2`;
- global Euclidean closest approach: `g*=0.5897545123`,
  `Delta_MC,min=0.5378414487`.

Finite-bandwidth protected-sector scan:

- 99 `(sigma_k/m, L)` points;
- 64 whole-class certificates;
- 33 explicit-overlap points;
- every explicit origin overlap is forced to the broad-class blind value
  `g_M=0`; this is not an overlap at the natural `g_M=2` vector point;
- 2 unresolved points: `(0.08,24)` and `(0.10,16)`;
- the positivity reduction plus pointwise charge envelope certifies 18 of the
  21 points unresolved by the earlier componentwise test;
- a 32-segment profile explicitly establishes overlap at `(0.15,8)`, where the
  eight-segment global QCQP minimum remains positive.
- the explicit Gaussian charge measure gives `<x>_C=0.1034492843` at
  `sigma_k/m=0.15`;
- the derived whole-class unit-covariance distance bound at `(0.10,2)` is
  `0.4042553595`;
- imposing the same numerical slope bound in Sachs variables gives a reference
  distance `0.4037104574` and 63 rather than 64 sufficient origin exclusions.

Matched detector-level fermion-versus-vector closure:

- fermion acceptance: `0.92758`;
- vector acceptance: `0.91766`;
- 12^3 + lost raw minimum KL: `0.1559913222`;
- estimated finite-N plug-in bias: `0.0070936675`;
- bias-corrected minimum KL: `0.1488976547`;
- recentered-bootstrap 95% interval: `[0.1447812963, 0.1529124751]`;
- classifier cross-check: `0.1736703 +/- 0.0040008`.

The collider calculation is a fixed-energy detector-readout closure test, not
a gauge-complete vector model and not a test of the analytic protected-sector
theorem. Its readout ablation is also not identified with the primitive
magnetic and charge-sensitive resources used in that theorem.

## Repository layout

```text
theory_numerics/
  finite_bandwidth_prd_v1_5.py
  detector_measures_distance_v1_5.py
  verify_publication_tables_v1_5.py
  results/prd_v1_5/   machine-readable v1.5 certificates and searches
  nuisance_domain_sensitivity.py
  legacy v1.1 verification and plotting scripts

collider/
  sufficient_statistics_v2.json.gz.b64
  scripts/              public reproducer and validation analyses
  results/              validated numerical outputs

generator_configuration/
  exact process definitions, PYTHIA command files, model metadata,
  software versions, detector-card provenance, and workflow provenance
```

## Reproduce the PRD v1.5 checkpoints

Python 3.12 is recommended.

```bash
python -m pip install -r requirements.txt

python theory_numerics/nuisance_domain_sensitivity.py
python theory_numerics/finite_bandwidth_prd_v1_5.py \
  --old-status theory_numerics/results/S1_slope_bandwidth_certificates_monotone_v1_1.csv \
  --outdir reproduced_results/prd_v1_5 \
  --unresolved-overlap-starts 1024

python theory_numerics/detector_measures_distance_v1_5.py \
  --outdir reproduced_results/prd_v1_5 \
  --figdir reproduced_results/figures

python theory_numerics/verify_publication_tables_v1_5.py

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

The GitHub Actions workflow runs the same public checks. It verifies the
explicit measure normalization, the derived distance and convergence tables,
the 64-versus-63 coordinate audit, the `64/33/2` classification and the two
unresolved-point search records, checks
the static nuisance-domain criterion, decodes the lossless collider sufficient
statistics, recomputes the public collider checkpoints, and syntax-checks the
event-level analysis scripts.

The finite-bandwidth classification is conservative. A positive analytic
lower bound proves separation over the full anchored Lipschitz class, and an
explicit zero proves overlap. Failed counterexample or certificate searches
are never promoted to proofs; this is why two scan points remain unresolved.

## Generator and detector provenance

The validated private production workflow used:

- MadGraph5_aMC@NLO 3.7.2 pinned commit
  `be1e7b273ca961c335ff2ee6da3688b5049b069e`;
- PYTHIA 8.312;
- ROOT 6.40.02;
- Delphes 3.5.1 with the CLICdet Stage-1 card;
- exact integer seeds and analysis commits in the machine-readable metadata.

The exact final process definitions and PYTHIA command files are in
`generator_configuration/`. The private `MadGraphPythiaDelphesDV01`
development repository is useful for production and CI, but is not required
to verify the publication-facing numerical claims.

## Citation

See `CITATION.cff` and cite the fixed GitHub Release `v1.5.0`. A separate
Zenodo deposit is not required for this repository's version-specific
reproducibility citation.
