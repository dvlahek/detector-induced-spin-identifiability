#!/usr/bin/env python3
"""Finite-sample validation for the reconstructed fermion-vs-vector KL diagnostic.

The primary estimator is the same discrete full detector channel used by the
closure analysis, with an explicit lost-event category. We calibrate the
finite-N plug-in bias by conditional multinomial resampling around the
high-statistics empirical reference, construct bootstrap intervals, scan
visible binning, and compare to an out-of-fold classifier likelihood-ratio
cross-check.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOG5 = math.log(5.0)


def kl(p: np.ndarray, q: np.ndarray) -> float:
    p = np.maximum(np.asarray(p, float), 1e-300)
    q = np.maximum(np.asarray(q, float), 1e-300)
    return float(np.sum(p * np.log(p / q)))


def directional(p: np.ndarray, q: np.ndarray) -> tuple[float, float, float]:
    fv = kl(p, q)
    vf = kl(q, p)
    return fv, vf, min(fv, vf)


def counts3(x: np.ndarray, bins: int) -> np.ndarray:
    edges = [np.linspace(0.0, 1.0, bins + 1)] * 3
    h, _ = np.histogramdd(x, bins=edges)
    return h.astype(np.int64).ravel()


def estimate_full(counts: np.ndarray, total: int, pseudo: float) -> np.ndarray:
    counts = np.asarray(counts, float)
    accepted = float(counts.sum())
    shape = (counts + pseudo) / (accepted + pseudo * len(counts))
    acceptance = accepted / float(total)
    out = np.concatenate([acceptance * shape, [1.0 - acceptance]])
    return out / out.sum()


def empirical_full(counts: np.ndarray, total: int) -> np.ndarray:
    lost = int(total - int(np.sum(counts)))
    if lost < 0:
        raise ValueError("Negative lost-event count")
    return np.concatenate([counts.astype(float), [float(lost)]]) / float(total)


def reestimate(sample: np.ndarray, pseudo: float) -> np.ndarray:
    sample = np.asarray(sample, dtype=np.int64)
    return estimate_full(sample[:-1], int(sample.sum()), pseudo)


def classifier_crosscheck(fermion_obs: np.ndarray, vector_obs: np.ndarray, Af: float, Av: float, seed: int, max_per_class: int) -> dict[str, object]:
    try:
        from sklearn.ensemble import HistGradientBoostingClassifier
        from sklearn.model_selection import StratifiedKFold
    except Exception as exc:
        return {"status": "skipped", "reason": str(exc)}

    rng = np.random.default_rng(seed)
    n = min(len(fermion_obs), len(vector_obs), max_per_class)
    fi = rng.choice(len(fermion_obs), size=n, replace=False)
    vi = rng.choice(len(vector_obs), size=n, replace=False)
    xf = np.asarray(fermion_obs[fi], dtype=np.float64)
    xv = np.asarray(vector_obs[vi], dtype=np.float64)
    X = np.vstack([xf, xv])
    y = np.concatenate([np.ones(n, dtype=np.int8), np.zeros(n, dtype=np.int8)])

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    folds = []
    for fold, (train, test) in enumerate(cv.split(X, y), start=1):
        clf = HistGradientBoostingClassifier(
            loss="log_loss", learning_rate=0.06, max_iter=180,
            max_leaf_nodes=31, min_samples_leaf=40,
            l2_regularization=1.0, early_stopping=True,
            random_state=seed + fold,
        )
        clf.fit(X[train], y[train])
        eta = np.clip(clf.predict_proba(X[test])[:, 1], 1e-5, 1.0 - 1e-5)
        logr = np.log(eta / (1.0 - eta))
        yt = y[test]
        dshape_fv = float(np.mean(logr[yt == 1]))
        dshape_vf = float(np.mean(-logr[yt == 0]))

        bern_fv = Af * math.log(Af / Av) + (1.0 - Af) * math.log((1.0 - Af) / (1.0 - Av))
        bern_vf = Av * math.log(Av / Af) + (1.0 - Av) * math.log((1.0 - Av) / (1.0 - Af))
        dfull_fv = Af * dshape_fv + bern_fv
        dfull_vf = Av * dshape_vf + bern_vf
        folds.append({
            "fold": fold,
            "shape_D_fermion_to_vector": dshape_fv,
            "shape_D_vector_to_fermion": dshape_vf,
            "full_D_fermion_to_vector": dfull_fv,
            "full_D_vector_to_fermion": dfull_vf,
            "full_D_min": min(dfull_fv, dfull_vf),
        })

    vals = np.asarray([x["full_D_min"] for x in folds], float)
    return {
        "status": "ok",
        "method": "5-fold out-of-fold HistGradientBoosting classifier odds with equal class priors",
        "events_per_class": n,
        "full_D_min_mean": float(np.mean(vals)),
        "full_D_min_sd_across_folds": float(np.std(vals, ddof=1)),
        "folds": folds,
        "interpretation": "Approximate binning-independent likelihood-ratio cross-check; not the primary estimator.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pseudo", type=float, default=0.5)
    ap.add_argument("--binnings", type=int, nargs="+", default=[6, 8, 10, 12])
    ap.add_argument("--replicates", type=int, default=300)
    ap.add_argument("--bootstrap-replicates", type=int, default=500)
    ap.add_argument("--classifier-max-per-class", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=20260831)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = np.load(args.reference)
    fobs = np.asarray(data["fermion_obs"], dtype=np.float64)
    vobs = np.asarray(data["vector_obs"], dtype=np.float64)
    nf = int(data["fermion_generated"])
    nv = int(data["vector_generated"])
    if nf != nv:
        raise RuntimeError("Validation expects equal generated sample sizes")
    nref = nf
    Af = len(fobs) / nf
    Av = len(vobs) / nv
    rng = np.random.default_rng(args.seed)

    sample_sizes = sorted(set([n for n in [10000, 20000, 50000, 100000, 300000, 500000, 1000000] if n <= nref] + [nref]))
    rows = []
    bias_rows = []

    for bins in args.binnings:
        cf = counts3(fobs, bins)
        cv = counts3(vobs, bins)
        pref_f = estimate_full(cf, nf, args.pseudo)
        pref_v = estimate_full(cv, nv, args.pseudo)
        reference = np.asarray(directional(pref_f, pref_v), dtype=float)
        bias_at_ref = np.full(3, np.nan)

        for n in sample_sizes:
            draws = np.empty((args.replicates, 3), dtype=float)
            for r in range(args.replicates):
                sf = rng.multinomial(n, pref_f)
                sv = rng.multinomial(n, pref_v)
                draws[r] = directional(reestimate(sf, args.pseudo), reestimate(sv, args.pseudo))
            bias = np.mean(draws, axis=0) - reference
            if n == nref:
                bias_at_ref = bias
            for j, direction in enumerate(["fermion_to_vector", "vector_to_fermion", "min"]):
                bias_rows.append({
                    "bins_per_axis": bins,
                    "sample_size_per_hypothesis": n,
                    "direction": direction,
                    "reference_D": reference[j],
                    "mean_estimated_D": float(np.mean(draws[:, j])),
                    "bias": float(bias[j]),
                    "sd": float(np.std(draws[:, j], ddof=1)),
                })

        if np.any(~np.isfinite(bias_at_ref)):
            raise RuntimeError("Could not calibrate bias at reference sample size")

        emp_f = empirical_full(cf, nf)
        emp_v = empirical_full(cv, nv)
        boot = np.empty((args.bootstrap_replicates, 3), dtype=float)
        for r in range(args.bootstrap_replicates):
            sf = rng.multinomial(nf, emp_f)
            sv = rng.multinomial(nv, emp_v)
            boot[r] = directional(reestimate(sf, args.pseudo), reestimate(sv, args.pseudo))

        corrected = reference - bias_at_ref
        center = np.mean(boot, axis=0)
        lo = corrected + (np.quantile(boot, 0.025, axis=0) - center)
        hi = corrected + (np.quantile(boot, 0.975, axis=0) - center)
        dmin = float(corrected[2])
        dmin_hi = float(hi[2])
        rows.append({
            "bins_per_axis": bins,
            "visible_categories": bins**3,
            "reference_events_per_hypothesis": nref,
            "fermion_acceptance": Af,
            "vector_acceptance": Av,
            "raw_full_D_fermion_to_vector": reference[0],
            "raw_full_D_vector_to_fermion": reference[1],
            "raw_full_D_min": reference[2],
            "estimated_bias_at_reference_min": bias_at_ref[2],
            "bias_corrected_full_D_fermion_to_vector": corrected[0],
            "bias_corrected_full_D_vector_to_fermion": corrected[1],
            "bias_corrected_full_D_min": corrected[2],
            "bias_corrected_D_min_ci025": lo[2],
            "bias_corrected_D_min_ci975": hi[2],
            "BH_necessary_events_point": LOG5 / dmin if dmin > 0 else math.inf,
            "BH_necessary_events_conservative": LOG5 / dmin_hi if dmin_hi > 0 else math.inf,
        })

    bias_df = pd.DataFrame(bias_rows)
    bin_df = pd.DataFrame(rows)
    bias_df.to_csv(args.outdir / "fermion_vector_kl_bias_calibration.csv", index=False)
    bin_df.to_csv(args.outdir / "fermion_vector_kl_binning_stability.csv", index=False)

    target_bins = 12 if 12 in args.binnings else args.binnings[-1]
    target = bin_df[bin_df["bins_per_axis"] == target_bins].iloc[0].to_dict()
    clf = classifier_crosscheck(fobs, vobs, Af, Av, args.seed + 77, args.classifier_max_per_class)
    summary = {
        "benchmark": "matched fermion-vs-W-like-vector antler closure",
        "reference_events_per_hypothesis": nref,
        "target_bins_per_axis": target_bins,
        "target": target,
        "classifier_crosscheck": clf,
        "interpretation": "Primary result is the bias-corrected discrete detector KL. The classifier is an approximate independent cross-check. BH values are necessary lower bounds only.",
    }
    (args.outdir / "fermion_vector_kl_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    q = bias_df[(bias_df["bins_per_axis"] == target_bins) & (bias_df["direction"] == "min")].copy()
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.plot(1.0 / q["sample_size_per_hypothesis"], q["bias"], marker="o")
    ax.set_xlabel("1 / events per hypothesis")
    ax.set_ylabel("estimated plug-in bias in min KL")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.outdir / "fermion_vector_kl_bias_vs_inverse_N.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.plot(bin_df["bins_per_axis"], bin_df["raw_full_D_min"], marker="o", label="raw")
    ax.plot(bin_df["bins_per_axis"], bin_df["bias_corrected_full_D_min"], marker="s", label="bias corrected")
    ax.set_xlabel("visible bins per axis")
    ax.set_ylabel("min reconstructed KL")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.outdir / "fermion_vector_kl_binning_stability.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
