# v1.4.0 — PRD reproducibility snapshot

This fixed release accompanies the PRD-focused manuscript **Whole-Class
Certificates for Detector-Restricted Identifiability of Elementary-Particle
Spin**.

## Central numerical update

- 99 anchored-Lipschitz scan points.
- 64 whole-class certificates.
- 33 explicit overlaps.
- 2 unresolved points, retained as unresolved rather than inferred from a
  failed search.
- The positivity reduction and pointwise charge envelope certify 18 of the 21
  points left unresolved by the earlier componentwise test.
- A 32-segment profile gives an explicit overlap at
  `(sigma_k/m, L) = (0.15, 8)`.
- Independent 32-segment searches with 1024 random starts remain positive at
  `(0.08, 24)` and `(0.10, 16)`; these are search records, not certificates.

## Reproducibility contents

The release source archive contains the complete public code and processed
data needed to reproduce the publication-facing numerical checkpoints:

- whole-Lipschitz, finite-basis QCQP, SOS, overlap, and unresolved-search
  records;
- static nuisance-domain sensitivity calculations;
- collider sufficient statistics and validated readout tables;
- exact MadGraph/PYTHIA process definitions and software provenance.

The private `MadGraphPythiaDelphesDV01` repository remains a development and CI
workspace. It is not required for the publication-facing reproduction path.

## Citation

Please cite this release as version `v1.4.0` using the repository
`CITATION.cff`. The permanent version-specific URL is:

https://github.com/dvlahek/detector-induced-spin-identifiability/releases/tag/v1.4.0
