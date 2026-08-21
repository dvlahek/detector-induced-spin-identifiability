# Finite-bandwidth protected-sector validation

The PRD v1.5 calculation uses the dimensionless transfer variable
`x = Q^2/m^2`, static conditions `G1(0)=1`, `GM(0)=g_M`,
`GQ(0)=1-g_M`, and a common anchored Lipschitz bound `|dGi/dx| <= L`.

The high-resolution grid contains 99 points:

- `sigma_k/m = 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50`;
- `L = 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24`.

The fixed primary classification contains 64 whole-class certificates, 33
explicit overlaps, and 2 unresolved points. The positivity reduction implies
that an origin overlap requires `G_M` to vanish. The remaining charge sector is
then bounded below by a pointwise envelope over the full anchored Lipschitz
class. This certifies 18 of the 21 points left unresolved in v1.1. Magnetic
positivity also forces every constructive origin overlap to `g_M=0`.

`finite_bandwidth_prd_v1_5.py` also performs three independent finite-basis
audits:

- eight-segment global QCQP minima by vertex enumeration and convex box QPs;
- directly verifiable degree-2 SOS identities;
- 32-segment searches, including the explicit overlap at `(0.15,8)` and the
  two unresolved searches at `(0.08,24)` and `(0.10,16)`.

These finite-basis computations do not replace the whole-class argument. A
failed search is not a proof, and the pointwise envelope is a conservative
relaxation because its minimizer need not form one globally compatible
Lipschitz function.

The companion `detector_measures_distance_v1_5.py` defines the Gaussian
detector measures explicitly, derives the whole-class covariance-whitened
half-strip enclosure, evaluates quadrature convergence, and repeats the
sufficient test in Sachs coordinates. At `(sigma_k/m,L)=(0.10,2)`, the primary
unit-covariance distance bound is `0.4042553595`; the Sachs-coordinate audit is
`0.4037104574`. The exact-origin envelope certifies 64 primary points and 63
Sachs-coordinate points.

Run:

```bash
python theory_numerics/finite_bandwidth_prd_v1_5.py \
  --old-status theory_numerics/results/S1_slope_bandwidth_certificates_monotone_v1_1.csv \
  --outdir reproduced_results/prd_v1_5 \
  --unresolved-overlap-starts 1024

python theory_numerics/detector_measures_distance_v1_5.py \
  --outdir reproduced_results/prd_v1_5 \
  --figdir reproduced_results/figures

python theory_numerics/verify_publication_tables_v1_5.py
```

The archived outputs are in `results/prd_v1_5/`. The older scripts and tables
remain available as historical inputs and independent cross-checks.
