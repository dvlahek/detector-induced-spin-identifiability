#!/usr/bin/env python3
"""Reproduce the publication estimator-space recovery audit.

The cloud-averaged columns use the explicitly stated nuisance distribution:
g is uniform on [-4,4], and each of the 24 slopes in the three eight-segment
form factors is independently uniform on [-2,2].  They are not minimax power
bounds.  A separate column evaluates the fixed admissible profile that supplies
the upper endpoint of the reference class-distance bracket.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


SEED = 20260816
SIGMA = 0.10
SLOPE_BOUND = 2.0
N_SEGMENTS = 8
GRID_N = 2601
KMAX_SIGMA = 9.0
INFERENCE_CLOUD_SIZE = 30000
TRUTH_CLOUD_SIZE = 10000
MC_DRAWS = 8000
EPSILON_LEVELS = (0.50, 0.30, 0.20, 0.15)
NEAR_BOUNDARY_FINGERPRINT = np.array([
    0.2677438808009649,
    0.30289269495740834,
])


def spectral_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return k, x, and the normalized magnetic and charge measures."""
    k = np.linspace(-KMAX_SIGMA * SIGMA, KMAX_SIGMA * SIGMA, GRID_N)
    energy = np.sqrt(1.0 + k * k)
    x = 2.0 * (energy - 1.0)
    eta = x / 4.0
    rho0 = np.exp(-0.5 * (k / SIGMA) ** 2) / (2.0 * energy)
    rho_m = rho0 * k * k
    rho_c = rho0 * (k / (energy + 1.0)) ** 2 * eta
    nu_m = rho_m / np.trapezoid(rho_m, k)
    nu_c = rho_c / np.trapezoid(rho_c, k)
    return k, x, nu_m, nu_c


def piecewise_profile(x: np.ndarray, slopes: np.ndarray,
                      anchor: float) -> np.ndarray:
    xmax = float(np.max(x))
    knots = np.linspace(0.0, xmax, N_SEGMENTS + 1)
    delta = np.diff(knots)
    starts = np.empty(len(knots))
    starts[0] = anchor
    starts[1:] = anchor + np.cumsum(slopes * delta)
    index = np.searchsorted(knots, x, side="right") - 1
    index = np.clip(index, 0, N_SEGMENTS - 1)
    return starts[index] + slopes[index] * (x - knots[index])


def fingerprint(parameters: np.ndarray, spectral: tuple[np.ndarray, ...]) -> np.ndarray:
    k, x, nu_m, nu_c = spectral
    g = float(parameters[0])
    slopes_1 = parameters[1:1 + N_SEGMENTS]
    slopes_m = parameters[1 + N_SEGMENTS:1 + 2 * N_SEGMENTS]
    slopes_q = parameters[1 + 2 * N_SEGMENTS:1 + 3 * N_SEGMENTS]
    g_1 = piecewise_profile(x, slopes_1, 1.0)
    g_m = piecewise_profile(x, slopes_m, g)
    g_q = piecewise_profile(x, slopes_q, 1.0 - g)
    w_m = np.trapezoid(nu_m * g_m * g_m, k)
    w_c = np.trapezoid(nu_c * (g_1 * g_q + 0.25 * x * g_q * g_q), k)
    return np.array([w_m, w_c], dtype=float)


def model_cloud(size: int, rng: np.random.Generator,
                spectral: tuple[np.ndarray, ...]) -> np.ndarray:
    cloud = np.empty((size, 2), dtype=float)
    for index in range(size):
        g = rng.uniform(-4.0, 4.0)
        slopes = rng.uniform(-SLOPE_BOUND, SLOPE_BOUND,
                             size=3 * N_SEGMENTS)
        cloud[index] = fingerprint(np.concatenate([[g], slopes]), spectral)
    return cloud


def roc_auc(null_scores: np.ndarray, alternative_scores: np.ndarray) -> float:
    combined = np.concatenate([null_scores, alternative_scores])
    ranks = pd.Series(combined).rank(method="average").to_numpy()
    n0 = len(null_scores)
    n1 = len(alternative_scores)
    return float((ranks[n0:].sum() - n1 * (n1 + 1) / 2) / (n0 * n1))


def run() -> pd.DataFrame:
    spectral = spectral_grid()
    rng = np.random.default_rng(SEED + 3333)
    inference = model_cloud(INFERENCE_CLOUD_SIZE, rng, spectral)
    tree = cKDTree(inference)
    truth = model_cloud(TRUTH_CLOUD_SIZE, rng, spectral)
    near_boundary_rng = np.random.default_rng(SEED + 4444)

    rows: list[dict[str, float]] = []
    for epsilon in EPSILON_LEVELS:
        null_scores = np.empty(MC_DRAWS)
        cloud_scores = np.empty(MC_DRAWS)
        for draw in range(MC_DRAWS):
            observed_null = rng.multivariate_normal(
                np.zeros(2), epsilon * epsilon * np.eye(2))
            distance_null = tree.query(observed_null, k=1)[0]
            null_scores[draw] = (
                observed_null @ observed_null - distance_null * distance_null)

            true_fingerprint = truth[rng.integers(0, len(truth))]
            observed_cloud = true_fingerprint + rng.multivariate_normal(
                np.zeros(2), epsilon * epsilon * np.eye(2))
            distance_cloud = tree.query(observed_cloud, k=1)[0]
            cloud_scores[draw] = (
                observed_cloud @ observed_cloud - distance_cloud * distance_cloud)

        threshold = float(np.quantile(null_scores, 0.95))
        near_boundary_observations = (
            NEAR_BOUNDARY_FINGERPRINT
            + near_boundary_rng.multivariate_normal(
                np.zeros(2), epsilon * epsilon * np.eye(2), size=MC_DRAWS)
        )
        origin_distance_sq = np.einsum(
            "ij,ij->i", near_boundary_observations, near_boundary_observations)
        model_distance = tree.query(near_boundary_observations, k=1)[0]
        near_boundary_scores = origin_distance_sq - model_distance * model_distance

        rows.append({
            "epsilon": epsilon,
            "cloud_averaged_spin1_power": float(np.mean(cloud_scores > threshold)),
            "near_boundary_profile_power": float(
                np.mean(near_boundary_scores > threshold)),
            "cloud_averaged_roc_auc": roc_auc(null_scores, cloud_scores),
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("theory_numerics/results/S3_precision_based_recovery.csv"),
    )
    args = parser.parse_args()
    frame = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
