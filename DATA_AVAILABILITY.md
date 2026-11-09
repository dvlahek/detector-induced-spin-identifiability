# Code and data availability

The code and processed numerical results corresponding to the manuscript are available at the immutable [scientific snapshot](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/93c70ed482443730c1abc6c2891732fb1bd5bd4f).

The primary calculations are implemented in `theory_numerics/spin_alpha0_*.py`. Reference datasets under [`theory_numerics/results/alpha0/`](theory_numerics/results/alpha0/) contain the 99-point classification, 33 constructed overlap paths, covariance-distance bounds, moment-deformation scan, QCQP/SOS results, estimator-space power, 48-segment searches, saved profiles and quadrature convergence.

The [reproduction guide](docs/REPRODUCIBILITY.md) gives commands for Windows, Linux and macOS. No private event-level collider data are required. The results are conditional on the explicitly specified detector model and do not establish a complete scattering apparatus.
