# Data availability

The public repository contains the code and processed numerical data used for the principal results in the manuscript, including the finite-bandwidth certificates, explicit overlap construction, distance bounds, compact-prior magnetic margins, estimator-space checks, and the sufficient statistics used for the reported collider histogram and readout-ablation results.

For the collider example, `collider/sufficient_statistics_v2.json.gz.b64` is a lossless archive of the count-level sufficient statistics used by the public reproducer. The repository also records the process definitions, software versions, seeds, detector provenance, and model-construction metadata under `generator_configuration/`.

Two intermediate event-level artifacts are not part of the public archive:

- the reconstructed event sample used for the out-of-fold classifier cross-check;
- the derived vector UFO used for the auxiliary hard-event generation chain.

The public reproducibility package was frozen at the processed sufficient-statistics level used for the reported collider histogram and readout results. The classifier and hard-event-generation checks are therefore documented as auxiliary validation records and are not used to establish the principal theoretical conclusions. Their scripts, metadata, seeds, and derived numerical summaries are public.

If the journal editors or peer reviewers require the non-public intermediate artifacts during review, they can be provided privately by the author without undue restriction or delay.
