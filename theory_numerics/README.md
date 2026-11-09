# Numerical implementation

The scripts in this directory support the spin-1/2 versus spin-1 detector-response analysis. The primary charge model uses the squared longitudinal-current measure $d\nu_C\propto\rho_0r^2x\,dr$ and $A_C=1$.

## Entry point

From the repository root:

```bash
python reproduce.py --mode smoke
python reproduce.py --mode full
python reproduce.py --mode highres
```

- **Smoke** checks the saved primary classification and independently integrates four saved 48-segment profiles.
- **Full** regenerates the primary finite-width tables, covariance audits, QCQP/SOS checks, power calculations and figures.
- **Highres** runs the optional seeded 48-segment optimization. Its positive results do not establish separation.

Reference data: [`results/alpha0/`](results/alpha0/). Generated outputs: `reproduced_alpha0/` (ignored by Git).

## Evidence and assumptions

The analytic magnetic and charge-envelope bounds cover the stated anchored Lipschitz class. The finite-basis QCQP, numerical SOS and multistart searches are complementary checks, not replacements for those analytic bounds. A negative response gives a constructive overlap witness after its sign is verified. A positive local optimum leaves separation undecided.

See the [evidence map](../docs/RESULTS.md) and [full reproduction procedure](../docs/REPRODUCIBILITY.md).
