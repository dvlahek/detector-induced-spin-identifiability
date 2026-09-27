# Electromagnetic response and finite-bandwidth calculations

This directory contains the active calculations for the protected spin-1/2 versus spin-1 response model. The dimensionless transfer is \(x=Q^2/m^2\). Form factors obey the static anchors \(G_1(0)=1\), \(G_M(0)=g_M\) and \(G_Q(0)=1-g_M\), together with the anchored Lipschitz condition \(|G_i(x)-G_i(y)|\le L|x-y|\).

## Compact prior: a magnetic-only certificate

For \(K=[1,3]\), the magnetic response has the lower bound

\[
W_M\ge\int \max(1-Lx,0)^2\,d\nu_M(x)>0.
\]

The strict positivity holds for every finite Gaussian bandwidth and finite \(L\); the numerical grid evaluates the size of this analytic margin. Run `narrow_prior_magnetic_audit.py` from the repository root. The checked 99-row output is `results/narrow_prior_magnetic_margins.csv`.

## Broad nuisance class

For \(g_M\in[-4,4]\), `finite_bandwidth.py` evaluates a whole-Lipschitz charge envelope and constructs explicit origin overlaps. The 99-point scan has 64 whole-class separation certificates, 33 explicit overlaps and two unresolved points. The global eight-segment QCQP, degree-2 SOS checks and 32-segment searches audit a finite basis; they are not substituted for whole-class proofs.

`detector_measures_distance.py` evaluates the normalized detector measures, the distance enclosure and the Sachs-coordinate check. `precision_recovery.py` provides the estimator-space audit. `verify_publication_tables.py` checks the tabulated publication results. Run the commands in the root README.

Active broad-class results are stored in `results/publication_v1_0/`. The input `results/S1_slope_bandwidth_certificates_monotone_v1_1.csv` is retained because the reproduction command uses it. Earlier scripts and numerical snapshots remain under `archive/development_snapshots/`.
