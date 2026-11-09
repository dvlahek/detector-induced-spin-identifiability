# Scientific results and corresponding evidence

The primary detector model has $A_C(k)=1$ and uses the charge measure $d\nu_C\propto\rho_0(r)r^2x\,dr$. See [Scientific scope](SCOPE.md) for physical assumptions and excluded claims.

| Scientific result | Reproducing code | Fixed data |
| --- | --- | --- |
| All-finite-bandwidth magnetic margin with a nonblind moment prior | [`narrow_prior_magnetic_audit.py`](../theory_numerics/narrow_prior_magnetic_audit.py) | [Magnetic margins](../theory_numerics/results/narrow_prior_magnetic_margins.csv) |
| **63** sufficient certificates, **33** overlaps, **3** unresolved | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [99-cell table](../theory_numerics/results/alpha0/classification_99_alpha0.csv) |
| 33 continuous-class signed overlap paths | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [Endpoint signs and roots](../theory_numerics/results/alpha0/explicit_overlaps_alpha0.csv) |
| Reference lower/upper distance bracket | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [Distance and normalization summary](../theory_numerics/results/alpha0/numerical_summary_alpha0.json) |
| **61** uniform certificates for $|\delta|\le0.10$ | [`spin_alpha0_compute.py`](../theory_numerics/spin_alpha0_compute.py) | [Deformation scan](../theory_numerics/results/alpha0/delta_alpha0.csv) |
| 21 finite-basis QCQPs and **12** positive numerical SOS lower bounds | [`spin_alpha0_qcqp.py`](../theory_numerics/spin_alpha0_qcqp.py) | [QCQP/SOS table](../theory_numerics/results/alpha0/qcqp_sos_alpha0.csv), [numerical multipliers](../theory_numerics/results/alpha0/sos_witness_alpha0.json) |
| 48-segment adversarial resolution check | [Search script](../theory_numerics/alpha0_48segment_audit.py), [independent quadrature](../theory_numerics/verify_48segment_quadrature.py) | [Saved profiles](../theory_numerics/results/alpha0/alpha0_48segments_profiles.json), [response table](../theory_numerics/results/alpha0/alpha0_48segments.csv), [quadrature](../theory_numerics/results/alpha0/alpha0_48segments_quadrature.csv) |
| Quadrature convergence | [Primary numeric routine](../theory_numerics/spin_alpha0_compute.py) | [Nine-point convergence table](../theory_numerics/results/alpha0/convergence_alpha0.csv) |
| Estimator-space power | [`spin_alpha0_power.py`](../theory_numerics/spin_alpha0_power.py) | [Estimator table](../theory_numerics/results/alpha0/estimator_power_alpha0.csv) |

**Interpretation.** Whole-class positivity results are sufficient certificates based on the analytic relaxation. Negative form-factor responses are constructive witnesses, with endpoint signs evaluated numerically. The finite-basis QCQP, SOS and multistart results are supplementary numerical audits. Positive multistart values do not settle the three unresolved cells. The SOS matrices have numerical PSD checks, not formal interval-arithmetic certification.

## Charge-acceptance sensitivity

The illustrative fixed-acceptance family $A_C^{(\alpha)}(k)=[2m/(E_k+m)]^{\alpha}$ gives these grid counts. The reference implementation uses $\alpha=0$.

| $\alpha$ | Certified | Constructed overlaps | Unresolved |
| ---: | ---: | ---: | ---: |
| **0** | **63** | **33** | **3** |
| 1 | 63 | 33 | 3 |
| 2 | 64 | 33 | 2 |
| 3 | 64 | 33 | 2 |

These four fixed choices do not imply a theorem uniform over arbitrary detector acceptances.
