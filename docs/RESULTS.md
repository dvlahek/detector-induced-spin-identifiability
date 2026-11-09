# Results and evidence (unfiltered charge model, v1.9.0)

| Finding | Reproducing script | Frozen data | Boundary |
| --- | --- | --- | --- |
| All-finite-bandwidth magnetic theorem under a nonblind moment prior | [`narrow_prior_magnetic_audit.py`](../theory_numerics/narrow_prior_magnetic_audit.py) | [`narrow_prior_magnetic_margins.csv`](../theory_numerics/results/narrow_prior_magnetic_margins.csv) | Analytic under stated positive measure and Lipschitz class |
| **63 certified / 33 overlap / 3 unresolved** over 99 cells ($\alpha=0$, $\delta=0$) | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [`classification_99_alpha0.csv`](../theory_numerics/results/alpha0/classification_99_alpha0.csv) | 63 sufficient whole-class bounds, 33 signed numerical witness paths; 3 neither |
| 33 two-kink constructed overlaps | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [`explicit_overlaps_alpha0.csv`](../theory_numerics/results/alpha0/explicit_overlaps_alpha0.csv) | Sign changes numerically evaluated; exact roots conditional on signs |
| Distance bracket $0.4034833471\le d_{\rm class}\le0.4034946602$ | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [`numerical_summary_alpha0.json`](../theory_numerics/results/alpha0/numerical_summary_alpha0.json) | Analytic half-strip enclosure plus explicit profile, numerical endpoint evaluation |
| $|\delta|\le0.10$: 61 certificates | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [`delta_alpha0.csv`](../theory_numerics/results/alpha0/delta_alpha0.csv) | Fixed primary measure; uniform over the stated nuisance interval |
| 21 finite-basis QCQPs; at least 12 positive numerical SOS checks | [`spin_alpha0_qcqp.py`](../theory_numerics/spin_alpha0_qcqp.py) | [`qcqp_sos_alpha0.csv`](../theory_numerics/results/alpha0/qcqp_sos_alpha0.csv), [`sos_witness_alpha0.json`](../theory_numerics/results/alpha0/sos_witness_alpha0.json) | Finite basis and floating-point PSD checks, not whole-class separation or interval proof |
| Estimator-space power | [`spin_alpha0_power.py`](../theory_numerics/spin_alpha0_power.py) | [`estimator_power_alpha0.csv`](../theory_numerics/results/alpha0/estimator_power_alpha0.csv) | Specified nuisance cloud, sampling law, reference covariance; not event-level detector predictions |

## Fixed charge-acceptance sensitivity

Only the primary $\alpha=0$ row is used in the article's headline results.

| $\alpha$ | Certified | Explicit overlaps | Unresolved |
| ---: | ---: | ---: | ---: |
| **0 (primary, unfiltered)** | **63** | **33** | **3** |
| 1 | 63 | 33 | 3 |
| 2 (former primary) | 64 | 33 | 2 |
| 3 | 64 | 33 | 2 |

These are four independent fixed detector architectures, not a uniform-in-$\alpha$ bound. Legacy filtered tables under `results/publication_v1_0/` are kept for archival comparisons only.
