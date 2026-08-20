# Frozen collider validation

The public collider package records the publication-facing sufficient statistics and validated outputs exported from the frozen MadGraph5_aMC@NLO -> PYTHIA 8.312 -> Delphes 3.5.1 CLICdet-Stage1 benchmark.

The three accepted-event observables are

1. `|cos(theta_mu-)|`,
2. `|cos(theta_mu+)|`,
3. `(E_mu- + E_mu+)/sqrt(s)`.

The primary full-channel estimator uses a 3D histogram plus one explicit lost-event category. The detector-readout degradation study conditions on the number of accepted events and therefore treats the overall production rate and accepted yield as nuisance parameters.

## Public sufficient statistics

`collider/sufficient_statistics_v2.json.gz.b64` is a lossless text-safe archive of the histogram counts required by the public reproducer. The GitHub Actions workflow decodes it automatically and checks the primary reconstructed KL, the detector-readout ablations, and the central-acceptance near-degeneracy checkpoint.

The full event-level analysis scripts are retained in `scripts/` and syntax-checked by CI. The event-level reconstructed sample itself is not needed to reproduce the publication-facing histogram statistics. The validated out-of-fold classifier result, which does require event-level observables, is preserved as an archived cross-check in `results/fermion_vector_kl_validation_summary.json`.

The exact generator and detector provenance is documented in `../generator_configuration/`. The private development repository is not part of this public publication package.
