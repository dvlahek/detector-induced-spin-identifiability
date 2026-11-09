# Research artifact changelog

Version fields indicate scientific code stages, not necessarily tagged GitHub releases. For quantitative reproduction cite the **frozen SHA**.

## 1.9.0 — 2026-11-09 (current J. Phys. A submission candidate)

- Switched the primary charge-channel model to the **unfiltered squared-current measure** $d\nu_C\propto\rho_0 r^2x\,dr$ ($A_C=1$, $\alpha=0$).
- Recalculated the main 99-point classification (63/33/3), $\delta$ robustness, class-distance bracket, covariance and Sachs audits, estimator-space power, 21 QCQP minima and numerical SOS feasibility checks.
- Rebuilt figures and manuscript source, removed IOP editorial placeholders, corrected bibliography to frozen alpha=0 commit.
- Preserved the old filtered-charge results without rewriting their history.
- Frozen scientific code/data: [`46271df7ba4b`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/f711505b03b1a01ad4d4be8a50ce4d46dbda2080).

## 1.8.0 — 2026-11-09 (archived $\alpha=2$ model)

The earlier filtered charge acceptance $A_C=[2m/(E_k+m)]^2$ gave 64/33/2 on the 99-point grid. Its scientific data were frozen at [`5c48258f`](https://github.com/dvlahek/detector-induced-spin-identifiability/tree/5c48258f69758818e3ed524cb4d7c1e9d5697811). They are not the primary results of the current manuscript.

## Earlier development

Releases tagged v1.0.0, v1.4.0, v1.5.0, and publication-v1.0.x preserve historical PRX/PRD-era materials and do not imply journal acceptance.
