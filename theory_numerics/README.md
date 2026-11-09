# Electromagnetic response and finite-bandwidth calculations

This directory contains the active calculations for the protected spin-1/2 versus spin-1 response model. The dimensionless transfer is \(x=Q^2/m^2\). The Haberzettl current class has static anchors \(G_1(0)=1\), \(G_M(0)=g_M\) and \(G_Q(0)=1-g_M\). The robustness extension writes \(g_M+g_Q=1+\delta\), so \(G_Q(0)=1+\delta-g_M\), with the same anchored Lipschitz condition \(|G_i(x)-G_i(y)|\le L|x-y|\).

## Compact prior: a magnetic-only certificate

For \(K=[1,3]\), the magnetic response has the lower bound

\[
W_M\ge\int \max(1-Lx,0)^2\,d\nu_M(x)>0.
\]

The strict positivity holds for every finite Gaussian bandwidth and finite \(L\), and it is independent of \(\delta\). The numerical grid evaluates the size of this analytic margin. Run `narrow_prior_magnetic_audit.py` from the repository root. The 99-row reference output is [`results/narrow_prior_magnetic_margins.csv`](results/narrow_prior_magnetic_margins.csv).

## Broad nuisance class

For \(g_M\in[-4,4]\), `finite_bandwidth.py` evaluates a whole-Lipschitz charge envelope and constructs explicit origin overlaps. The 99-point scan has 64 whole-class separation certificates, 33 explicit overlaps and two unresolved points. The global eight-segment QCQP, degree-2 SOS checks and 32-segment searches audit a finite basis; they are not substituted for whole-class proofs.

`delta_robustness_audit.py` enlarges the static relation to \(g_M+g_Q=1+\delta\) and computes whole-class bounds uniform over \(|\delta|\le\Delta\). The summary is stored in [`results/delta_robustness_summary.csv`](results/delta_robustness_summary.csv). `detector_measures_distance.py` evaluates the normalized detector measures, the distance enclosure and the Sachs-coordinate check. `precision_recovery.py` provides the estimator-space audit. `verify_publication_tables.py` checks the tabulated publication results. Run the commands in the root README.

The broad-class reference data are stored in `results/publication_v1_0/`. This directory name identifies the fixed numerical dataset used by the manuscript; it is not a separate, superseded manuscript. The input `results/S1_slope_bandwidth_certificates_monotone_v1_1.csv` is kept alongside the active scripts because the finite-bandwidth reproduction uses its earlier status labels. Superseded scripts and earlier result sets are under [`../archive/development_snapshots/`](../archive/development_snapshots/).


## Version 1.8.0 additional audits

A chosen positive charge acceptance `A_C(k)=(2m/(E_k+m))**2` explains the normalized measure used by the manuscript. Run `charge_kernel_identity_check.py` to check the measure identity. Run `explicit_overlap_family_verify.py` to reproduce the signed endpoint witness for all 33 overlap cells. Run `acceptance_sensitivity_audit.py` to calculate four separate acceptance choices alpha=0,1,2,3; counts (certified, overlap, unresolved) are (63,33,3), (63,33,3), (64,33,2), (64,33,2). These sensitivity checks do not prove uniformity over an interval of alpha.
