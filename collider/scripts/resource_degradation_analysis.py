#!/usr/bin/env python3
"""Detector-readout degradation study for the frozen 100k fermion-vector closure.

The analysis uses only the accepted reconstructed observables stored in the
validated 100k artifact. It conditions on a fixed number N of accepted
opposite-sign dimuon events, so total production rate and acceptance yield are
nuisance parameters. Nested readout architectures are compared by a
Jeffreys-smoothed discrete KL diagnostic, conditional finite-N bias calibration,
bootstrap uncertainty, equal-prior likelihood-ratio pseudoexperiments, and an
uncertainty-aware Pinsker lower bound on the optimal binary error.

The readout ablation is a collider-level detector-architecture study. It is not
identified with the primitive magnetic/charge resources in the analytic
protected-sector theorem.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def kl(p, q):
    p = np.clip(np.asarray(p, dtype=float), 1e-300, None)
    q = np.clip(np.asarray(q, dtype=float), 1e-300, None)
    return float(np.sum(p * np.log(p / q)))


def directional(p, q):
    a, b = kl(p, q), kl(q, p)
    return np.asarray([a, b, min(a, b)], dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--bins", type=int, default=12)
    ap.add_argument("--pseudo", type=float, default=0.5)
    ap.add_argument("--bias-reps", type=int, default=500)
    ap.add_argument("--bootstrap-reps", type=int, default=1000)
    ap.add_argument("--mc-trials", type=int, default=50000)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    dat = np.load(args.reference)
    f = np.asarray(dat["fermion_obs"], dtype=float)
    v = np.asarray(dat["vector_obs"], dtype=float)
    edges = np.linspace(0.0, 1.0, args.bins + 1)

    architectures = {
        "full_spectrometric_3D": (0, 1, 2),
        "direction_pair_2D": (0, 1),
        "single_direction_1D": (0,),
        "energy_only_1D": (2,),
    }
    labels = {
        "full_spectrometric_3D": "angles + energy",
        "direction_pair_2D": "two directions",
        "single_direction_1D": "one direction",
        "energy_only_1D": "energy only",
    }

    def counts(arr, dims, eta_max):
        zc = math.tanh(eta_max)
        x = arr[(arr[:, 0] < zc) & (arr[:, 1] < zc)]
        if len(dims) == 1:
            h = np.histogram(x[:, dims[0]], bins=edges)[0]
        else:
            h = np.histogramdd(x[:, dims], bins=[edges] * len(dims))[0].ravel()
        return h.astype(np.int64), int(len(x))

    def prob(c):
        c = np.asarray(c, dtype=float)
        return (c + args.pseudo) / (c.sum() + args.pseudo * len(c))

    def calibrate(dims, eta_max, seed, bias_reps=None, boot_reps=None):
        br = args.bias_reps if bias_reps is None else bias_reps
        rr = args.bootstrap_reps if boot_reps is None else boot_reps
        cf, nf = counts(f, dims, eta_max)
        cv, nv = counts(v, dims, eta_max)
        pf, pv = prob(cf), prob(cv)
        raw = directional(pf, pv)
        rng = np.random.default_rng(seed)

        draws = np.empty((br, 3), dtype=float)
        for i in range(br):
            draws[i] = directional(
                prob(rng.multinomial(nf, pf)),
                prob(rng.multinomial(nv, pv)),
            )
        bias = draws.mean(axis=0) - raw
        corrected = raw - bias

        ef, ev = cf / cf.sum(), cv / cv.sum()
        boot = np.empty((rr, 3), dtype=float)
        for i in range(rr):
            boot[i] = directional(
                prob(rng.multinomial(nf, ef)),
                prob(rng.multinomial(nv, ev)),
            )
        center = boot.mean(axis=0)
        lo = corrected + np.quantile(boot, 0.025, axis=0) - center
        hi = corrected + np.quantile(boot, 0.975, axis=0) - center
        return {
            "nf": nf,
            "nv": nv,
            "pf": pf,
            "pv": pv,
            "raw": raw,
            "bias": bias,
            "corrected": corrected,
            "lo": lo,
            "hi": hi,
        }

    def lr_error(p, q, n, trials, seed):
        """Bayes error for the empirical discrete model, estimated by MC."""
        rng = np.random.default_rng(seed)
        llr = np.log(np.clip(p, 1e-300, None) / np.clip(q, 1e-300, None))
        errs = []
        for dist, positive in ((p, True), (q, False)):
            wrong = 0.0
            done = 0
            while done < trials:
                m = min(1000, trials - done)
                idx = rng.choice(len(dist), size=(m, n), replace=True, p=dist)
                score = llr[idx].sum(axis=1)
                if positive:
                    wrong += np.sum(score < 0) + 0.5 * np.sum(score == 0)
                else:
                    wrong += np.sum(score > 0) + 0.5 * np.sum(score == 0)
                done += m
            errs.append(wrong / trials)
        pe = 0.5 * (errs[0] + errs[1])
        se = 0.5 * math.sqrt(
            errs[0] * (1 - errs[0]) / trials
            + errs[1] * (1 - errs[1]) / trials
        )
        return pe, se, errs

    def pinsker_lower(d_upper, n):
        """Pe* >= 1/2(1-sqrt(N D/2)); use an upper limit for D."""
        return max(
            0.0,
            0.5 * (1.0 - math.sqrt(max(0.0, n * max(0.0, d_upper) / 2.0))),
        )

    rows, calibrated = [], {}
    for i, (name, dims) in enumerate(architectures.items()):
        c = calibrate(dims, 2.5, args.seed + 1000 * i)
        calibrated[name] = c
        rows.append({
            "architecture": name,
            "label": labels[name],
            "eta_max": 2.5,
            "fermion_events": c["nf"],
            "vector_events": c["nv"],
            "raw_D_min": c["raw"][2],
            "estimated_bias_min": c["bias"][2],
            "bias_corrected_D_min": c["corrected"][2],
            "ci95_low": c["lo"][2],
            "ci95_high": c["hi"][2],
        })
    arch = pd.DataFrame(rows)
    arch.to_csv(args.outdir / "architecture_kl_summary.csv", index=False)

    budgets = [1, 5, 10, 20, 50, 100, 200]
    err_rows = []
    for ai, (name, c) in enumerate(calibrated.items()):
        d95_upper = float(max(0.0, c["hi"][2]))
        for n in budgets:
            pe, se, errs = lr_error(
                c["pf"], c["pv"], n, args.mc_trials,
                args.seed + 100000 * ai + n,
            )
            err_rows.append({
                "architecture": name,
                "label": labels[name],
                "accepted_event_budget_N": n,
                "plug_in_Bayes_error": pe,
                "MC_standard_error": se,
                "fermion_error": errs[0],
                "vector_error": errs[1],
                "bias_corrected_D95_upper": d95_upper,
                "Pinsker_lower_bound_from_D95_upper": pinsker_lower(d95_upper, n),
            })
    errors = pd.DataFrame(err_rows)
    errors.to_csv(args.outdir / "bayes_error_vs_event_budget.csv", index=False)

    eta_values = [0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5]
    eta_rows = []
    for j, eta in enumerate(eta_values):
        c = calibrate((0, 1), eta, args.seed + 50000 + 1000 * j)
        d95_upper = float(max(0.0, c["hi"][2]))
        for n in budgets:
            eta_rows.append({
                "eta_max": eta,
                "accepted_event_budget_N": n,
                "fermion_events_available": c["nf"],
                "vector_events_available": c["nv"],
                "raw_D_min": c["raw"][2],
                "bias_corrected_D_min": c["corrected"][2],
                "D_min_ci95_low": c["lo"][2],
                "D_min_ci95_high": c["hi"][2],
                "Pinsker_lower_bound_from_D95_upper": pinsker_lower(d95_upper, n),
            })
    eta_df = pd.DataFrame(eta_rows)
    eta_df.to_csv(args.outdir / "direction_only_acceptance_resource_map.csv", index=False)

    barrel = calibrate(
        (0, 1), 0.8, args.seed + 8800,
        bias_reps=2000, boot_reps=5000,
    )
    barrel_d95_upper = float(max(0.0, barrel["hi"][2]))
    barrel_summary = {
        "architecture": "direction_pair_2D",
        "eta_max": 0.8,
        "fermion_events_in_acceptance": barrel["nf"],
        "vector_events_in_acceptance": barrel["nv"],
        "raw_D_min": float(barrel["raw"][2]),
        "estimated_bias_min": float(barrel["bias"][2]),
        "bias_corrected_D_min_unconstrained": float(barrel["corrected"][2]),
        "physical_point_estimate_after_nonnegative_boundary": float(max(0.0, barrel["corrected"][2])),
        "recentered_bootstrap_ci95_low_unconstrained": float(barrel["lo"][2]),
        "recentered_bootstrap_ci95_high": float(barrel["hi"][2]),
        "bias_calibration_replicates": 2000,
        "bootstrap_replicates": 5000,
    }
    (args.outdir / "barrel_direction_highstat_verification.json").write_text(
        json.dumps(barrel_summary, indent=2), encoding="utf-8"
    )

    barrel_err = []
    for n in budgets:
        pe, se, _ = lr_error(
            barrel["pf"], barrel["pv"], n, args.mc_trials,
            args.seed + 880000 + n,
        )
        barrel_err.append({
            "accepted_event_budget_N": n,
            "plug_in_Bayes_error": pe,
            "MC_standard_error": se,
            "bias_corrected_D95_upper": barrel_d95_upper,
            "Pinsker_lower_bound_from_D95_upper": pinsker_lower(barrel_d95_upper, n),
        })
    pd.DataFrame(barrel_err).to_csv(
        args.outdir / "barrel_direction_bayes_error.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(7.3, 4.9))
    for name in architectures:
        sub = errors[errors.architecture == name]
        ax.errorbar(
            sub.accepted_event_budget_N,
            sub.plug_in_Bayes_error,
            yerr=sub.MC_standard_error,
            marker="o",
            capsize=2,
            label=labels[name],
        )
    ax.axhline(0.5, linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_ylim(0, 0.52)
    ax.set_xlabel("accepted event budget $N$")
    ax.set_ylabel("plug-in equal-prior Bayes error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.outdir / "bayes_error_vs_event_budget.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    pivot = errors.pivot(
        index="architecture",
        columns="accepted_event_budget_N",
        values="plug_in_Bayes_error",
    )
    order = list(architectures)
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    im = ax.imshow(pivot.loc[order, budgets].to_numpy(), aspect="auto", origin="upper")
    ax.set_xticks(range(len(budgets)), labels=[str(x) for x in budgets])
    ax.set_yticks(range(len(order)), labels=[labels[x] for x in order])
    ax.set_xlabel("accepted event budget $N$")
    fig.colorbar(im, ax=ax, label="plug-in $P_e^*$")
    fig.tight_layout()
    fig.savefig(args.outdir / "resource_bayes_error_map.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    piv2 = eta_df.pivot(
        index="eta_max",
        columns="accepted_event_budget_N",
        values="Pinsker_lower_bound_from_D95_upper",
    )
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    im = ax.imshow(piv2.loc[eta_values, budgets].to_numpy(), aspect="auto", origin="lower")
    ax.set_xticks(range(len(budgets)), labels=[str(x) for x in budgets])
    ax.set_yticks(range(len(eta_values)), labels=[str(x) for x in eta_values])
    ax.set_xlabel("accepted event budget $N$")
    ax.set_ylabel(r"direction-only acceptance $|\eta|<\eta_{\max}$")
    fig.colorbar(im, ax=ax, label=r"lower bound on $P_e^\star$")
    fig.tight_layout()
    fig.savefig(args.outdir / "direction_only_conservative_error_map.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "study_scope": "collider-level readout ablation; not the primitive magnetic/charge theorem resources",
        "conditioning": "shape-only at fixed accepted event count N; total rate and accepted yield are nuisance parameters",
        "baseline_acceptance": "pT_mu>10 GeV, |eta_mu|<2.5",
        "bins_per_retained_coordinate": args.bins,
        "jeffreys_pseudocount": args.pseudo,
        "bias_replicates": args.bias_reps,
        "bootstrap_replicates": args.bootstrap_reps,
        "likelihood_ratio_trials_per_hypothesis": args.mc_trials,
        "pinsker_convention": "uses upper 95% recentered bootstrap limit of bias-corrected minimum directional KL",
        "seed": args.seed,
        "barrel_highstat": barrel_summary,
    }
    (args.outdir / "resource_study_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(arch.to_string(index=False))
    print(json.dumps(barrel_summary, indent=2))


if __name__ == "__main__":
    main()
