#!/usr/bin/env python3
"""Reproduce public collider checkpoints from histogram sufficient statistics.

The event-level MadGraph/PYTHIA/Delphes development chain is intentionally not
part of this public repository. The archived histogram counts are sufficient
to reproduce the primary discrete KL diagnostic, readout-ablation raw KLs, and
independent finite-sample bias/bootstrap calibrations.

The out-of-fold classifier cross-check requires event-level observables and is
therefore archived as a validated output rather than recomputed here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def kl(p, q):
    p = np.maximum(np.asarray(p, float), 1e-300)
    q = np.maximum(np.asarray(q, float), 1e-300)
    return float(np.sum(p * np.log(p / q)))


def directional(p, q):
    a = kl(p, q)
    b = kl(q, p)
    return np.asarray([a, b, min(a, b)], dtype=float)


def estimate_full(counts, total, pseudo=0.5):
    counts = np.asarray(counts, float)
    accepted = float(counts.sum())
    shape = (counts + pseudo) / (accepted + pseudo * len(counts))
    acceptance = accepted / float(total)
    out = np.concatenate([acceptance * shape, [1.0 - acceptance]])
    return out / out.sum()


def shape_prob(counts, pseudo=0.5):
    counts = np.asarray(counts, float)
    return (counts + pseudo) / (counts.sum() + pseudo * len(counts))


def reestimate_full(sample, pseudo=0.5):
    sample = np.asarray(sample, np.int64)
    return estimate_full(sample[:-1], int(sample.sum()), pseudo)


def fresh_full_calibration(cf, cv, nf, nv, seed, reps=300, boot_reps=500):
    pf = estimate_full(cf, nf)
    pv = estimate_full(cv, nv)
    raw = directional(pf, pv)
    rng = np.random.default_rng(seed)
    draws = np.empty((reps, 3), float)
    for r in range(reps):
        sf = rng.multinomial(nf, pf)
        sv = rng.multinomial(nv, pv)
        draws[r] = directional(reestimate_full(sf), reestimate_full(sv))
    bias = draws.mean(axis=0) - raw
    corrected = raw - bias

    lost_f = nf - int(np.sum(cf))
    lost_v = nv - int(np.sum(cv))
    emp_f = np.concatenate([np.asarray(cf, float), [lost_f]]) / nf
    emp_v = np.concatenate([np.asarray(cv, float), [lost_v]]) / nv
    boot = np.empty((boot_reps, 3), float)
    for r in range(boot_reps):
        sf = rng.multinomial(nf, emp_f)
        sv = rng.multinomial(nv, emp_v)
        boot[r] = directional(reestimate_full(sf), reestimate_full(sv))
    center = boot.mean(axis=0)
    lo = corrected + np.quantile(boot, 0.025, axis=0) - center
    hi = corrected + np.quantile(boot, 0.975, axis=0) - center
    return raw, bias, corrected, lo, hi


def fresh_shape_calibration(cf, cv, nf, nv, seed, reps=500, boot_reps=1000):
    pf, pv = shape_prob(cf), shape_prob(cv)
    raw = directional(pf, pv)
    rng = np.random.default_rng(seed)
    draws = np.empty((reps, 3), float)
    for r in range(reps):
        draws[r] = directional(
            shape_prob(rng.multinomial(nf, pf)),
            shape_prob(rng.multinomial(nv, pv)),
        )
    bias = draws.mean(axis=0) - raw
    corrected = raw - bias
    ef = np.asarray(cf, float) / np.sum(cf)
    ev = np.asarray(cv, float) / np.sum(cv)
    boot = np.empty((boot_reps, 3), float)
    for r in range(boot_reps):
        boot[r] = directional(
            shape_prob(rng.multinomial(nf, ef)),
            shape_prob(rng.multinomial(nv, ev)),
        )
    center = boot.mean(axis=0)
    lo = corrected + np.quantile(boot, 0.025, axis=0) - center
    hi = corrected + np.quantile(boot, 0.975, axis=0) - center
    return raw, bias, corrected, lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--statistics", type=Path, required=True)
    ap.add_argument("--full12", type=Path, required=True)
    ap.add_argument("--archived-closure", type=Path, required=True)
    ap.add_argument("--archived-resource", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    d = json.loads(args.statistics.read_text())
    full = json.loads(args.full12.read_text())
    nf = int(d["generated"]["fermion"])
    nv = int(d["generated"]["vector"])

    cf = np.asarray(full["fermion_counts"], np.int64)
    cv = np.asarray(full["vector_counts"], np.int64)
    raw, bias, corr, lo, hi = fresh_full_calibration(
        cf, cv, nf, nv, seed=20260825, reps=300, boot_reps=500
    )

    assert abs(raw[2] - 0.1559913221919423) < 1e-12
    assert 0.005 < bias[2] < 0.0095
    assert 0.145 < corr[2] < 0.153
    assert lo[2] > 0.135

    expected_raw = {
        "full_spectrometric_3D": 0.16922142194653927,
        "direction_pair_2D": 0.0062287072843407315,
        "single_direction_1D": 0.0012738276442316315,
        "energy_only_1D": 0.152866426933505,
    }
    resource = {}
    for i, name in enumerate(expected_raw):
        h = full if name == "full_spectrometric_3D" else d["histograms"][name]
        cfi = np.asarray(h["fermion_counts"], np.int64)
        cvi = np.asarray(h["vector_counts"], np.int64)
        nfi = int(h["fermion_accepted"])
        nvi = int(h["vector_accepted"])
        r, b, c, l, u = fresh_shape_calibration(
            cfi, cvi, nfi, nvi, seed=20260817 + 1000 * i
        )
        assert abs(r[2] - expected_raw[name]) < 1e-12
        resource[name] = {
            "raw_D_min": float(r[2]),
            "fresh_bias_min": float(b[2]),
            "fresh_corrected_D_min": float(c[2]),
            "fresh_ci95": [float(l[2]), float(u[2])],
        }

    assert 0.0035 < resource["direction_pair_2D"]["fresh_corrected_D_min"] < 0.0060
    assert 0.158 < resource["full_spectrometric_3D"]["fresh_corrected_D_min"] < 0.164

    barrel = d["direction_eta_scan"]["0.8"]
    br, bb, bc, bl, bu = fresh_shape_calibration(
        np.asarray(barrel["fermion_counts"], np.int64),
        np.asarray(barrel["vector_counts"], np.int64),
        int(barrel["fermion_accepted"]),
        int(barrel["vector_accepted"]),
        seed=20260817 + 8800,
        reps=2000,
        boot_reps=5000,
    )
    assert abs(br[2] - 0.002288681067189543) < 1e-12
    assert abs(bc[2]) < 2.5e-4
    assert bu[2] < 0.003

    archived_closure = json.loads(args.archived_closure.read_text())
    target = archived_closure["target"]
    assert abs(target["bias_corrected_full_D_min"] - 0.1488976546638089) < 1e-15
    assert abs(target["bias_corrected_D_min_ci025"] - 0.1447812962754047) < 1e-15
    assert abs(target["bias_corrected_D_min_ci975"] - 0.15291247510327222) < 1e-15

    archived_resource = pd.read_csv(args.archived_resource).set_index("architecture")
    assert abs(archived_resource.loc["direction_pair_2D", "bias_corrected_D_min"] - 0.004677054428747011) < 1e-15
    assert abs(archived_resource.loc["full_spectrometric_3D", "bias_corrected_D_min"] - 0.16128791447730945) < 1e-15

    out = {
        "primary_full_channel": {
            "raw_D_min": float(raw[2]),
            "fresh_bias_min": float(bias[2]),
            "fresh_corrected_D_min": float(corr[2]),
            "fresh_ci95": [float(lo[2]), float(hi[2])],
            "archived_corrected_D_min": target["bias_corrected_full_D_min"],
            "archived_ci95": [target["bias_corrected_D_min_ci025"], target["bias_corrected_D_min_ci975"]],
        },
        "resource_architectures": resource,
        "barrel_direction_eta_0p8": {
            "raw_D_min": float(br[2]),
            "fresh_bias_min": float(bb[2]),
            "fresh_corrected_D_min": float(bc[2]),
            "fresh_ci95": [float(bl[2]), float(bu[2])],
        },
        "note": "Fresh Monte Carlo calibration uses public sufficient statistics. Exact archived publication values are independently checked against frozen output tables.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2))
    print("public_collider_statistics_ok")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
