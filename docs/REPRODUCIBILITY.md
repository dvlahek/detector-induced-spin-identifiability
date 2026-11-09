# Reproduce the numerical results

## Environment and immutable source

Use Python 3.12 and the [frozen scientific commit](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/2421638e1311dd2319a269d5c98dcfa94cb57970) when reproducing the fixed numerical results. `requirements.txt` specifies dependency ranges rather than a complete platform lockfile. Retain the output of `python -m pip freeze` when reporting a cross-platform discrepancy.

```bash
git clone https://github.com/dvlahek/detector-induced-spin-identifiability.git
cd detector-induced-spin-identifiability
git checkout 2421638e1311dd2319a269d5c98dcfa94cb57970
python -m pip install -r requirements.txt
```

## Run the cited snapshot

```bash
python theory_numerics/spin_alpha0_compute.py
python theory_numerics/spin_alpha0_qcqp.py
python theory_numerics/spin_alpha0_power.py
python theory_numerics/spin_alpha0_figures.py
python theory_numerics/verify_alpha0_release.py
python theory_numerics/verify_48segment_quadrature.py
```

The commands regenerate the primary 99-point classification, seven $\Delta$ scans, covariance bounds, QCQP/SOS audits, estimator-space power and figures, then verify the saved tables and independently integrate the 48-segment witness profiles. On the current default branch, the convenience wrapper `python reproduce.py --mode full` performs the same primary workflow.

The scripts write to `theory_numerics/reproduced_alpha0/` and keep the repository reference data unchanged. The regenerated figures use Matplotlib. Power calculations use a fixed pseudorandom seed and specified finite inference and truth clouds.

## Optional high-resolution search

```bash
python theory_numerics/alpha0_48segment_audit.py --out theory_numerics/reproduced_alpha0 --starts 200
```

This uses a seeded 48-segment alternating optimization with 200 starts for the known overlap $(0.15,8)$ and the three unresolved cells $(0.08,24)$, $(0.10,16)$ and $(0.50,1)$. The frozen profile arrays and independent quadrature are stored in `theory_numerics/results/alpha0/alpha0_48segments_profiles.json` and `alpha0_48segments_quadrature.csv`. It may find a negative witness at the known-overlap control, but positive optimized values do not establish separation.

## Verification criteria

- 99 grid cells consist of 63 whole-class sufficient certificates, 33 signed overlap constructions and 3 unresolved cases.
- At $|\delta|\le0.10$, 61 grid cells retain positive whole-class bounds.
- The reference distance is enclosed between 0.4034833471 and 0.4034946602.
- The numerical degree-2 SOS audit has 12 positive saved feasible lower bounds among 21 finite-basis QCQPs.
- The 48-segment negative control is negative when independently integrated. The three unresolved cases are not promoted to positive theorems.

The SOS checks use floating-point feasibility, not interval arithmetic.

