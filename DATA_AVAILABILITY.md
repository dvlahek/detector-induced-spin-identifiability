# Code and numerical data

The complete code and processed data for the spin-identifiability calculations are publicly available at the [frozen scientific source](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/e4bb51f9e6891c88b46b8583b3fd0ef5c6414f8b).

The reference implementation is in [`theory_numerics/`](theory_numerics/). The primary datasets under [`theory_numerics/results/alpha0/`](theory_numerics/results/alpha0/) include 99-point classification, 33 constructive overlap paths, covariance-distance bounds, moment-deformation scans, QCQP/SOS results, estimator-space power, 48-segment profiles and quadrature convergence.

Run the independent verification or complete reproduction using the [reproducibility guide](docs/REPRODUCIBILITY.md). Generated data are written to a separate ignored directory. No private event-level collider data are required for the primary results.

The conclusions are conditional on the specified detector-response model and do not establish a complete experimental scattering apparatus.
