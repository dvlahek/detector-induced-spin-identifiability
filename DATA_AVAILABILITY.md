# Data and code availability

The present J. Phys. A submission candidate uses the **unfiltered longitudinal charge response** ($\alpha=0$). The reproducibility scripts and processed numerical data are available at the fixed commit [`46271df7ba4b94cfa6ad2af7f2627209c1e0ae7b`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/46271df7ba4b94cfa6ad2af7f2627209c1e0ae7b).

The primary evidence is under `theory_numerics/results/alpha0/` with calculation scripts `theory_numerics/spin_alpha0_*.py`. It includes 99-point classification, 33 signed witnesses, moment-deformation scan, covariance-distance bounds, QCQP/SOS numerical checks, and estimator-space power.

The earlier `theory_numerics/results/publication_v1_0/` and old scripts reproduce an **archival filtered-charge ($\alpha=2$)** model and are **not** the source for current manuscript values. Older code is retained for transparency and comparison.

All results relevant to this conditional theoretical manuscript are computed from public scripts and parameters; no non-public collider event records are required. Auxiliary collider folders are not a validation claim for the current manuscript.

No formally interval-certified numerical sign or completed scattering apparatus is claimed.
