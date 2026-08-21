# PRD reproducibility snapshot v1.5.0

This fixed release closes the main reproducibility and theorem-scope issues
identified in the pre-submission review.

- Defines the Gaussian detector weight and both normalized transfer measures.
- Adds a derived whole-class covariance-whitened distance enclosure; the
  reference unit-covariance bound is `0.4042553595`.
- Adds quadrature and Gaussian-tail convergence records.
- Repeats the sufficient origin-exclusion test in Sachs variables: 64 primary
  covariant certificates versus 63 Sachs-coordinate certificates.
- States explicitly that all 33 constructive origin overlaps occur at the
  broad-class blind value `g_M=0`.
- Retains the conservative `64/33/2` primary classification and the two
  unresolved 32-segment search records.
- Adds a fail-fast v1.5 publication-table verifier.

The collider outputs remain the same frozen fixed-energy detector-readout
closure. They are not a gauge-complete electroweak prediction and are not used
to prove the protected-sector theorem.
