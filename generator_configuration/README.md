# Generator and detector provenance

The publication-facing repository deliberately does **not** expose the private development/infrastructure repository. The frozen reconstructed reference in `../collider/frozen_data/reco_reference_100k.npz` is the public analysis input.

The validated event-production chain used:

- MadGraph5_aMC@NLO 3.7.2 pinned to commit `be1e7b273ca961c335ff2ee6da3688b5049b069e`;
- PYTHIA 8.312;
- ROOT 6.40.02;
- Delphes 3.5.1;
- Delphes card `cards/delphes_card_CLICdet_Stage1.tcl` from the Delphes 3.5.1 release (Git blob SHA `97bdd3aee86f07c895b4b4c684840cecfdff4c11`).

Benchmark:

- sqrt(s) = 500 GeV;
- charged-parent mass = 200 GeV;
- invisible mass = 100 GeV;
- 100000 generated events per hypothesis;
- generator seed = 20260824;
- validation seed = 20260825;
- photon s-channel production only; Z and neutrino t-channel exchange excluded.

The fermion hypothesis uses `MSSM_SLHA2`. The vector hypothesis is a W-like benchmark constructed from a copied pinned Standard Model UFO, keeping the Standard Model gamma-W-W Yang-Mills vertex while promoting the muon-neutrino field to a massive stable invisible state. Auxiliary values MZ=300 GeV and GF=7.544120778283774e-7 GeV^-2 were used so the internal tree-level MW equals the matched 200 GeV parent mass. This is a benchmark construction, not a Standard Model cross-section prediction.

The exact process definitions, PYTHIA command files, model metadata, software versions, and workflow provenance are archived in this directory.
