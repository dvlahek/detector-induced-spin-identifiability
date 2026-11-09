# Main results and reproducibility evidence

This evidence map connects the current **J. Phys. A manuscript in preparation** to public code, data, and the precise scope of each claim.

| Result | Public source | Boundary |
| --- | --- | --- |
| Magnetic separation at every finite Gaussian bandwidth and finite $L$ with $\inf_{g\in K}|g|>0$ | [`narrow_prior_magnetic_audit.py`](../theory_numerics/narrow_prior_magnetic_audit.py), [99-row output](../theory_numerics/results/narrow_prior_magnetic_margins.csv), Proposition 7.2 | Requires a nonblind moment domain; not a universal physical prior |
| Broad reference grid: **64 certified, 33 constructive overlaps, 2 unresolved** | [`finite_bandwidth_status.csv`](../theory_numerics/results/publication_v1_0/finite_bandwidth_status.csv), [publication table verifier](../theory_numerics/verify_publication_tables.py) | Chosen charge measure, anchored Lipschitz class, $\delta=0$, $g_M\in[-4,4]$ |
| Constructive two-kink witnesses in all 33 overlap cells | [`explicit_overlap_family_verify.py`](../theory_numerics/explicit_overlap_family_verify.py), [signed witnesses](../theory_numerics/results/publication_v1_0/explicit_overlap_family_results.csv) | Endpoint signs numerically checked; exact root via continuity and IVT |
| Unresolved cells $(0.08,24)$, $(0.10,16)$ | [`unresolved_32segment_search.json`](../theory_numerics/results/publication_v1_0/unresolved_32segment_search.json) | Positive finite-basis minima cannot certify the entire function class |
| Reference unit-covariance class-distance bracket $0.4042553595\le d_{\rm class}\le0.4042657175$ | [`reference_distance_bracket.json`](../theory_numerics/results/publication_v1_0/reference_distance_bracket.json) and [`detector_measures_distance.py`](../theory_numerics/detector_measures_distance.py) | Analytic bound with numerically evaluated endpoints; no formal interval arithmetic |
| $|\delta|\le0.10$: 61 sufficient certificates | [`delta_robustness_audit.py`](../theory_numerics/delta_robustness_audit.py), [summary](../theory_numerics/results/delta_robustness_summary.csv) | Uniform over stated nuisance interval, not an experimental coupling bound |
| Chosen charge acceptance and four fixed alternatives | [`charge_kernel_identity_check.py`](../theory_numerics/charge_kernel_identity_check.py), [`acceptance_sensitivity_audit.py`](../theory_numerics/acceptance_sensitivity_audit.py) | Fixed-detector-model choices, not all possible acceptances |

## Charge-acceptance sensitivity

For $A_C^{(\alpha)}=[2m/(E_k+m)]^\alpha$:

| $\alpha$ | Certified | Overlap | Unresolved |
| ---: | ---: | ---: | ---: |
| 0 | 63 | 33 | 3 |
| 1 | 63 | 33 | 3 |
| 2 (reference) | 64 | 33 | 2 |
| 3 | 64 | 33 | 2 |

Complete fixed-choice data: [99-point by four CSV](../theory_numerics/results/acceptance_sensitivity_99point.csv). These four evaluations **do not establish uniform robustness for a continuous interval of $\alpha$**.

## Data locations

- [`theory_numerics/results/publication_v1_0/`](../theory_numerics/results/publication_v1_0/) retains publication-oriented result snapshots and QCQP/SOS audits. The name is **historical**, not a superseded data claim.
- [`theory_numerics/results/`](../theory_numerics/results/) contains additional v1.8.0, $\delta$ and compact-prior results.
- [`manuscript_jphysa/`](../manuscript_jphysa/) contains the article proof statements, assumptions, main/supplement source and figures.

Automated checks reproduce specified mathematical/numerical outputs. They are **not** an independent experimental replication or a demonstrated instrument.
