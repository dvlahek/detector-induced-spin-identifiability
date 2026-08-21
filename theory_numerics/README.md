# Publication 1.0 finite-bandwidth validation

The active calculation uses `x=Q^2/m^2`, static conditions
`G1(0)=1`, `GM(0)=g_M`, `GQ(0)=1-g_M`, and the anchored Lipschitz bound
`|dGi/dx| <= L`.

The 99-point grid contains nine values of `sigma_k/m` and eleven values of
`L`. Its fixed classification is 64 whole-class certificates, 33 explicit
overlaps, and 2 unresolved points. The calculation also supplies:

- normalized detector measures and their convergence audits;
- a deterministic whole-class half-strip distance enclosure on all 99 points;
- a tight lower/upper distance bracket at `(0.10,2)`;
- covariant-versus-Sachs sufficient-envelope comparison;
- a direct Sachs-slope audit of the explicit overlap witness;
- eight-segment QCQP and degree-2 SOS finite-basis audits;
- 1024-start, 32-segment searches at the two unresolved points.
- a nuisance-averaged estimator-space recovery audit plus a fixed
  near-boundary-profile power column.

Run the three commands in the repository README. Active archived outputs are
in `results/publication_v1_0/`. Superseded scripts and versioned results are in
`../archive/development_snapshots/theory_numerics/`.
