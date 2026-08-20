# Finite-bandwidth protected-sector validation

The manuscript studies finite detector bandwidth using the dimensionless transfer variable `x = Q^2/m^2`, static conditions `G1(0)=1`, `GM(0)=g_M`, `GQ(0)=1-g_M`, and a common Lipschitz bound `|dGi/dx| <= L`.

The high-resolution publication grid is:

- `g_M in [-4,4]`;
- `sigma_k/m = 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50`;
- `L = 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24`;
- 99 scan points;
- 8 piecewise-linear adversarial segments;
- 2601 spectral grid points;
- maximum sampled `|k| = 9 sigma_k`;
- differential-evolution maximum iterations 180;
- population multiplier 16;
- numerical zero tolerance `2e-4`;
- 8000 estimator Monte Carlo draws per hypothesis;
- model cloud size 30000.

The archived final status table contains 46 guaranteed points, 32 explicit-overlap points after monotone closure, and 21 unresolved points. The three-way classification is deliberately conservative: failure of the analytic certificate does not establish overlap, and failure of the adversarial search does not establish separation.

`verify_publication_tables.py` checks all publication-facing archived checkpoints. `plot_status_map.py` recreates the final status map. `monotone_closure.py` documents the post-processing rule for nested slope classes.
