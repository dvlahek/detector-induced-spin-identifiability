# v1.1.0 - nuisance-domain sensitivity revision

This revision adds the static nuisance-domain sensitivity analysis used in the revised manuscript.

## New publication checkpoint

For the elementary static protected fingerprint

`(W_M,W_C)=(g^2,1-g)`,

the minimum primitive protected support depends explicitly on the admitted nuisance class `K`:

- `tau(K)=1` if `0` is excluded from `K` or `1` is excluded from `K`;
- `tau(K)=2` only when both single-channel blind points `0` and `1` are admitted.

For the centered family `K_R=[2-R,2+R]`, the charge-only blind point enters at `R=1`, while the magnetic blind point enters at `R=2`. Therefore the two-resource requirement begins exactly at `R=2`. The combined Euclidean margin remains positive, with global minimum `0.5378414487` at `g=0.5897545123`.

## Added files

- `theory_numerics/nuisance_domain_sensitivity.py`
- `theory_numerics/results/nuisance_domain_sensitivity.csv`
- generated `theory_numerics/figures/nuisance_domain_sensitivity.png`

The GitHub Actions reproduction workflow now verifies these checkpoints in addition to the existing finite-bandwidth, collider, and detector-readout checks.
