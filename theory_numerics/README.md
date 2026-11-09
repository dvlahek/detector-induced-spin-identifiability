# Active theoretical and numerical calculations

## Current J. Phys. A analysis: **unfiltered $\alpha=0$**

The charge measure is $d\nu_C\propto\rho_0(r)r^2 x\,dr$ (normalized). A positive detector acceptance $A_C=1$ means the squared longitudinal-current response has **no additional momentum-dependent event filtering**. The magnetic measure remains unchanged. Analysis uses real anchored Lipschitz form factors and the specified $g_M+g_Q=1+\delta$ moment relation.

**Frozen code snapshot:** [`46271df7ba4b94cfa6ad2af7f2627209c1e0ae7b`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/46271df7ba4b94cfa6ad2af7f2627209c1e0ae7b).

```bash
python theory_numerics/spin_alpha0_compute.py
python theory_numerics/spin_alpha0_qcqp.py
python theory_numerics/spin_alpha0_power.py
python theory_numerics/spin_alpha0_figures.py
python theory_numerics/verify_alpha0_release.py
```

Current results: [`results/alpha0/`](results/alpha0/). Scripts regenerate in ignored `reproduced_alpha0/` and are checked against the committed snapshot. The 99 cells are **63 certified / 33 overlap / 3 unresolved**; at $|\delta|\le0.10$, 61 are certified. Numerically evaluated positivity and SOS matrices do **not** constitute interval-arithmetic certification.

## Previous filtered-charge stage (archive)

Files such as `finite_bandwidth.py`, `detector_measures_distance.py`, `precision_recovery.py`, `delta_robustness_audit.py` and `results/publication_v1_0/` refer to the **older $\alpha=2$ model**. They are retained for provenance only and **must not** be used as the current unfiltered-model evidence. The magnetic-only audit is acceptance-independent and remains applicable.

See [scientific scope](../docs/SCOPE.md) and [full reproducibility instructions](../docs/REPRODUCIBILITY.md).
