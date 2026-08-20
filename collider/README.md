# Frozen collider validation

The public collider analysis is based on `frozen_data/reco_reference_100k.npz`, exported from the validated full MadGraph5_aMC@NLO -> PYTHIA 8.312 -> Delphes 3.5.1 CLICdet-Stage1 chain.

The three accepted-event observables are

1. `|cos(theta_mu-)|`,
2. `|cos(theta_mu+)|`,
3. `(E_mu- + E_mu+)/sqrt(s)`.

The full primary estimator uses a 3D histogram plus one explicit lost-event category. The readout-degradation study conditions on the number of accepted events and therefore treats the overall cross section and accepted yield as nuisance parameters.

The scripts in `scripts/` are copied from the validated publication commit. The archived outputs in `results/` and `figures/` come from the successful GitHub Actions certificates used in the manuscript.
