# Generator and detector provenance

The public reproduction path uses
`../collider/sufficient_statistics_v2.json.gz.b64`. This losslessly encoded
file contains the histogram counts needed for every public collider checkpoint.
The event-level reconstructed sample and the modified vector UFO used during
development are not distributed here. Scripts that require those inputs are
retained as audit code, and their validated outputs are archived under
`../collider/results/`.

The validated event-production chain used:

- MadGraph5_aMC@NLO 3.7.2 pinned to commit `be1e7b273ca961c335ff2ee6da3688b5049b069e`;
- PYTHIA 8.312;
- ROOT 6.40.02;
- Delphes 3.5.1;
- the CLICdet Stage-1 card distributed as
  `cards/delphes_card_CLICdet_Stage1.tcl` in Delphes 3.5.1 (upstream Git blob
  SHA `97bdd3aee86f07c895b4b4c684840cecfdff4c11`; the upstream card is identified
  exactly but is not duplicated in this repository).

Benchmark:

- sqrt(s) = 500 GeV;
- charged-parent mass = 200 GeV;
- invisible mass = 100 GeV;
- 100000 generated events per hypothesis;
- fermion shower seed = 20260824;
- vector shower seed = 20260925;
- public sufficient-statistics validation seed = 20260825;
- event-level validation default seed = 20260831;
- detector-resource analysis seed = 20260817;
- photon s-channel production only; Z and neutrino t-channel exchange excluded.

The fermion hypothesis uses `MSSM_SLHA2`. The vector hypothesis is a W-like benchmark constructed from a copied pinned Standard Model UFO, keeping the Standard Model gamma-W-W Yang-Mills vertex while promoting the muon-neutrino field to a massive stable invisible state. Auxiliary values MZ=300 GeV and GF=7.544120778283774e-7 GeV^-2 were used so the internal tree-level MW equals the matched 200 GeV parent mass. This is a benchmark construction, not a Standard Model cross-section prediction.

This directory records the exact process definitions, PYTHIA command files,
model metadata, software versions, and workflow provenance used to produce the
published summaries. Re-running event generation additionally requires the
modified vector UFO and event-level sample noted above.
