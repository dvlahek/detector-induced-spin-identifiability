#!/usr/bin/env python3
"""PRD v1.5 detector-measure, distance, and coordinate-robustness audit.

This driver makes the two normalized detector measures explicit, evaluates a
whole-class rectangular enclosure of the protected fingerprint for every
static anchor g, minimizes its covariance-whitened distance from the origin,
and repeats the exact-origin charge-envelope test in covariant (G1,GQ) and
Sachs (GC,GQ) Lipschitz coordinates.

The reported distance is a lower bound because the pointwise extrema at
different x need not assemble into a single anchored Lipschitz profile.  The
formula is analytic; the numerical tables include grid/tail convergence so
that the quadrature evaluation can be audited independently.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


SIGMAS = (0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
SLOPES = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0)


@dataclass(frozen=True)
class DetectorMeasures:
    k: np.ndarray
    x: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray


def detector_measures(sigma: float, grid_n: int = 20001,
                      kmax_sigma: float = 12.0) -> DetectorMeasures:
    """Return the normalized measures for m=1 and a Gaussian spectral kernel.

    With E=sqrt(1+k^2), x=2(E-1), and
        r0=exp[-k^2/(2 sigma^2)]/(2E),
    the unnormalized densities in k are
        rho_M=r0 k^2,
        rho_C=r0 [k/(E+1)]^2 (x/4).
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    k = np.linspace(-kmax_sigma * sigma, kmax_sigma * sigma, grid_n)
    energy = np.sqrt(1.0 + k * k)
    x = 2.0 * (energy - 1.0)
    rho0 = np.exp(-0.5 * (k / sigma) ** 2) / (2.0 * energy)
    rho_m = rho0 * k * k
    rho_c = rho0 * (k / (energy + 1.0)) ** 2 * (x / 4.0)
    norm_m = np.trapezoid(rho_m, k)
    norm_c = np.trapezoid(rho_c, k)
    if norm_m <= 0.0 or norm_c <= 0.0:
        raise RuntimeError("non-positive detector normalization")
    return DetectorMeasures(k=k, x=x, nu_m=rho_m / norm_m,
                            nu_c=rho_c / norm_c)


def average(measures: DetectorMeasures, density: np.ndarray,
            values: np.ndarray) -> float:
    return float(np.trapezoid(density * values, measures.k))


def charge_pointwise_bounds(x: np.ndarray, g: float, slope: float,
                            quadratic_coefficient: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact rectangle extrema of a*b+c*x*b^2 at each x.

    a is centered at 1 and b at 1-g, with radius Lx.  The extrema occur at
    b endpoints, at the stationary point for either a endpoint, or at the
    b=0 kink where the minimizing a endpoint switches.
    """
    a_lo, a_hi = 1.0 - slope * x, 1.0 + slope * x
    b_lo, b_hi = 1.0 - g - slope * x, 1.0 - g + slope * x
    candidates: list[np.ndarray] = []
    for a in (a_lo, a_hi):
        candidates.append(a * b_lo + quadratic_coefficient * x * b_lo * b_lo)
        candidates.append(a * b_hi + quadratic_coefficient * x * b_hi * b_hi)
        stationary = np.zeros_like(x)
        nonzero = x > 1.0e-15
        np.divide(-a, 2.0 * quadratic_coefficient * x,
                  out=stationary, where=nonzero)
        stationary = np.clip(stationary, b_lo, b_hi)
        candidates.append(a * stationary
                          + quadratic_coefficient * x * stationary * stationary)
    candidates.append(np.where((b_lo <= 0.0) & (b_hi >= 0.0), 0.0, np.nan))
    stack = np.asarray(candidates)
    return np.nanmin(stack, axis=0), np.nanmax(stack, axis=0)


def fingerprint_rectangle(measures: DetectorMeasures, g: float, slope: float,
                          coordinate: str = "G1") -> tuple[float, float, float]:
    """Return M_L(g), C_L^-(g), C_L^+(g) for the whole Lipschitz class."""
    magnetic_floor = average(
        measures, measures.nu_m,
        np.maximum(abs(g) - slope * measures.x, 0.0) ** 2,
    )
    coefficient = 0.25 if coordinate == "G1" else 1.0 / 12.0
    charge_lo, charge_hi = charge_pointwise_bounds(
        measures.x, g, slope, coefficient
    )
    return (
        magnetic_floor,
        average(measures, measures.nu_c, charge_lo),
        average(measures, measures.nu_c, charge_hi),
    )


def covariance_shape(ratio: float, correlation: float) -> np.ndarray:
    if ratio <= 0.0 or abs(correlation) >= 1.0:
        raise ValueError("invalid covariance shape")
    return np.array([[1.0, correlation * ratio],
                     [correlation * ratio, ratio * ratio]])


def distance_to_half_strip(precision: np.ndarray, magnetic_floor: float,
                           charge_lo: float, charge_hi: float) -> float:
    """Minimize w^T A w for w_M>=M and C_lo<=w_C<=C_hi."""
    trial_c = [charge_lo, charge_hi]
    if charge_lo <= 0.0 <= charge_hi:
        trial_c.append(0.0)
    trial_c.append(float(np.clip(
        -precision[0, 1] * magnetic_floor / precision[1, 1],
        charge_lo, charge_hi,
    )))
    if abs(precision[0, 1]) > 1.0e-15:
        trial_c.append(float(np.clip(
            -magnetic_floor * precision[0, 0] / precision[0, 1],
            charge_lo, charge_hi,
        )))
    values = []
    for charge in trial_c:
        magnetic = max(
            magnetic_floor,
            -precision[0, 1] * charge / precision[0, 0],
        )
        point = np.array([magnetic, charge])
        values.append(float(point @ precision @ point))
    return float(np.sqrt(max(0.0, min(values))))


def class_distance(sigma: float, slope: float, ratio: float = 1.0,
                   correlation: float = 0.0, coordinate: str = "G1",
                   grid_n: int = 20001, kmax_sigma: float = 12.0) -> tuple[float, float]:
    measures = detector_measures(sigma, grid_n=grid_n,
                                 kmax_sigma=kmax_sigma)
    precision = np.linalg.inv(covariance_shape(ratio, correlation))

    def objective(g: float) -> float:
        rectangle = fingerprint_rectangle(measures, float(g), slope, coordinate)
        return distance_to_half_strip(precision, *rectangle)

    anchors = np.linspace(-4.0, 4.0, 401)
    values = np.asarray([objective(g) for g in anchors])
    best = int(np.argmin(values))
    lo = anchors[max(0, best - 2)]
    hi = anchors[min(len(anchors) - 1, best + 2)]
    result = minimize_scalar(objective, bounds=(lo, hi), method="bounded",
                             options={"xatol": 1.0e-11})
    return float(result.fun), float(result.x)


def origin_charge_envelope(sigma: float, slope: float,
                           coordinate: str) -> float:
    measures = detector_measures(sigma)
    lower, _ = charge_pointwise_bounds(
        measures.x, g=0.0, slope=slope,
        quadratic_coefficient=0.25 if coordinate == "G1" else 1.0 / 12.0,
    )
    return average(measures, measures.nu_c, lower)


def plot_coordinate_map(frame: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.5), sharey=True)
    for ax, coordinate, title in zip(
        axes, ("G1", "GC"),
        (r"Lipschitz in $(G_1,G_M,G_Q)$",
         r"Lipschitz in $(G_C,G_M,G_Q)$"),
    ):
        sub = frame[frame.coordinate == coordinate]
        pivot = sub.pivot(index="L", columns="sigma_k_over_m",
                          values="origin_envelope_positive").reindex(
                              index=SLOPES, columns=SIGMAS)
        ax.imshow(pivot.values.astype(int), origin="lower", aspect="auto",
                  cmap=plt.matplotlib.colors.ListedColormap(["#d9d9d9", "#58a66a"]),
                  vmin=0, vmax=1)
        ax.set_xticks(range(len(SIGMAS)), [f"{s:g}" for s in SIGMAS], rotation=45)
        ax.set_yticks(range(len(SLOPES)), [f"{ell:g}" for ell in SLOPES])
        ax.set_xlabel(r"$\sigma_k/m$")
        ax.set_title(title, fontsize=11)
    axes[0].set_ylabel(r"Lipschitz bound $L$")
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", markersize=9,
                   markerfacecolor="#58a66a", markeredgecolor="none",
                   label="origin excluded"),
        plt.Line2D([0], [0], marker="s", linestyle="", markersize=9,
                   markerfacecolor="#d9d9d9", markeredgecolor="none",
                   label="criterion inconclusive"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_static_response(output: Path) -> None:
    g = np.linspace(-0.5, 2.15, 600)
    wm, wc = g * g, 1.0 - g
    gstar = 0.5897545123
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.plot(wm, wc, lw=2.4, color="#2878b5",
            label=r"elementary spin 1: $(g^2,1-g)$")
    points = {
        "spin $1/2$\nrank-2 origin": (0.0, 0.0, "#d62728", (16, -42)),
        "$g=0$\nmagnetic blind point": (0.0, 1.0, "#2878b5", (16, 8)),
        "$g=1$\ncharge blind point": (1.0, 0.0, "#9467bd", (14, 12)),
        r"$g_\star=0.590$" "\nclosest approach": (gstar ** 2, 1.0 - gstar,
                                             "#ff8c1a", (30, 30)),
        "$g=2$\nYang--Mills point": (4.0, -1.0, "#2ca02c", (-118, 22)),
    }
    for label, (xv, yv, color, offset) in points.items():
        ax.scatter([xv], [yv], s=74, color=color, zorder=4)
        ax.annotate(label, (xv, yv), xytext=offset, textcoords="offset points",
                    fontsize=9.5, arrowprops={"arrowstyle": "-", "color": "0.35"})
    ax.axhline(0.0, color="0.75", lw=0.8)
    ax.axvline(0.0, color="0.75", lw=0.8)
    ax.set_xlim(-0.15, 4.55)
    ax.set_ylim(-1.25, 1.8)
    ax.set_xlabel(r"$W_M$")
    ax.set_ylabel(r"$W_C$")
    ax.legend(loc="upper right", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--figdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.figdir.mkdir(parents=True, exist_ok=True)

    diagnostics = []
    for sigma in SIGMAS:
        measures = detector_measures(sigma)
        diagnostics.append({
            "sigma_k_over_m": sigma,
            "mean_x_M": average(measures, measures.nu_m, measures.x),
            "mean_x_C": average(measures, measures.nu_c, measures.x),
            "mean_x_C_over_sigma_squared":
                average(measures, measures.nu_c, measures.x) / sigma ** 2,
            "normalization_M": average(measures, measures.nu_m,
                                        np.ones_like(measures.x)),
            "normalization_C": average(measures, measures.nu_c,
                                        np.ones_like(measures.x)),
        })
    pd.DataFrame(diagnostics).to_csv(
        args.outdir / "detector_measure_diagnostics_v1_5.csv", index=False)

    coordinates = []
    for coordinate in ("G1", "GC"):
        for sigma in SIGMAS:
            for slope in SLOPES:
                bound = origin_charge_envelope(sigma, slope, coordinate)
                coordinates.append({
                    "coordinate": coordinate,
                    "sigma_k_over_m": sigma,
                    "L": slope,
                    "origin_charge_lower_envelope": bound,
                    "origin_envelope_positive": bool(bound > 1.0e-10),
                })
    coordinate_frame = pd.DataFrame(coordinates)
    coordinate_frame.to_csv(
        args.outdir / "coordinate_robustness_v1_5.csv", index=False)
    plot_coordinate_map(
        coordinate_frame,
        args.figdir / "coordinate_convention_certificate_map.png",
    )

    covariance_rows = []
    for ratio in (0.5, 1.0, 2.0, 4.0):
        for correlation in (0.0, 0.5, -0.5):
            distance, anchor = class_distance(
                sigma=0.10, slope=2.0, ratio=ratio,
                correlation=correlation, coordinate="G1")
            covariance_rows.append({
                "sigmaC_over_sigmaM": ratio,
                "rho": correlation,
                "whole_class_whitened_lower_bound": distance,
                "minimizing_anchor_gM": anchor,
                "estimator_scale_for_3sigma": distance / 3.0,
                "estimator_scale_for_5sigma": distance / 5.0,
            })
    pd.DataFrame(covariance_rows).to_csv(
        args.outdir / "whole_class_distance_covariance_v1_5.csv", index=False)

    convergence = []
    for coordinate in ("G1", "GC"):
        for grid_n in (10001, 20001, 40001):
            for tail in (10.0, 12.0, 14.0):
                distance, anchor = class_distance(
                    sigma=0.10, slope=2.0, coordinate=coordinate,
                    grid_n=grid_n, kmax_sigma=tail)
                convergence.append({
                    "coordinate": coordinate,
                    "grid_n": grid_n,
                    "kmax_sigma": tail,
                    "whole_class_whitened_lower_bound": distance,
                    "minimizing_anchor_gM": anchor,
                })
    pd.DataFrame(convergence).to_csv(
        args.outdir / "distance_quadrature_convergence_v1_5.csv", index=False)

    plot_static_response(args.figdir / "response_space_static.png")

    summary = {
        "origin_envelope_positive_counts": {
            coordinate: int(coordinate_frame[
                coordinate_frame.coordinate == coordinate
            ].origin_envelope_positive.sum())
            for coordinate in ("G1", "GC")
        },
        "reference_G1_unit_covariance": covariance_rows[3],
        "reference_GC_unit_covariance": dict(zip(
            ("whole_class_whitened_lower_bound", "minimizing_anchor_gM"),
            class_distance(0.10, 2.0, coordinate="GC"),
        )),
        "sigma_0p15_mean_x_C": diagnostics[4]["mean_x_C"],
        "sigma_0p15_mean_x_C_over_sigma_squared":
            diagnostics[4]["mean_x_C_over_sigma_squared"],
    }
    (args.outdir / "detector_measures_distance_v1_5_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
