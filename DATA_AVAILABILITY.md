# Data availability

The current J. Phys. A manuscript uses the public scripts and processed numerical data under `theory_numerics/`.

These files reproduce the compact-prior magnetic margins, the broad-class finite-bandwidth classification, explicit overlap witnesses, distance enclosures, coordinate checks, estimator-space audits, and the robustness calculation for `g_M + g_Q = 1 + delta`.

No non-public event-level data are required for the results reported in the current manuscript.

The `collider/` and `generator_configuration/` directories are retained as auxiliary project provenance and are not part of the current manuscript's validation chain.


The v1.8.0 numerical snapshot includes `theory_numerics/explicit_overlap_family_verify.py`, `theory_numerics/acceptance_sensitivity_audit.py`, and `theory_numerics/charge_kernel_identity_check.py` together with their CSV/JSON outputs. The finite-bandwidth claims are conditional on the specified nonnegative charge-acceptance architecture; alternative alpha choices are shown separately as a sensitivity study. No experimental apparatus or collider detection claim is inferred from these checks.
