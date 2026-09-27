#!/usr/bin/env python3
"""Magnetic-only finite-bandwidth lower bound for the compact prior K=[1,3].

For |G_M(x)-g| <= L*x and g >= 1, the calibrated magnetic response obeys
W_M >= integral max(1-L*x,0)^2 dnu_M. The bound is positive for every
finite L and every finite Gaussian bandwidth sigma_k/m > 0.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from scipy.integrate import quad

SIGMAS = (0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
SLOPES = (0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24)
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "results" / "narrow_prior_magnetic_margins.csv"


def envelope(sigma: float, slope: float) -> float:
    """Evaluate the lower bound using the normalized Gaussian magnetic measure."""
    if sigma <= 0 or not math.isfinite(sigma) or slope < 0 or not math.isfinite(slope):
        raise ValueError("sigma must be positive and finite; L must be finite and nonnegative")
    if slope == 0:
        return 1.0
    kcut = math.sqrt((1 + 1 / (2 * slope)) ** 2 - 1)

    def density(k: float) -> float:
        return math.exp(-k * k / (2 * sigma * sigma)) * k * k / (2 * math.sqrt(1 + k * k))

    def integrand(k: float) -> float:
        x = 2 * (math.sqrt(1 + k * k) - 1)
        return density(k) * max(1 - slope * x, 0) ** 2

    denominator = quad(density, 0, max(12 * sigma, kcut),
                       epsabs=1e-13, epsrel=2e-13, limit=220)[0]
    numerator = quad(integrand, 0, kcut, epsabs=1e-14,
                     epsrel=2e-13, limit=220)[0]
    return numerator / denominator


def main(output: Path) -> None:
    rows = [(sigma, slope, envelope(sigma, slope))
            for sigma in SIGMAS for slope in SLOPES]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("sigma_k_over_m", "L", "K_low", "K_high",
                         "magnetic_only_lower_bound"))
        writer.writerows((s, L, 1, 3, format(value, ".14g")) for s, L, value in rows)
    minimum = min(rows, key=lambda row: row[2])
    print(f"K=[1,3]: {sum(v > 0 for _, _, v in rows)}/{len(rows)} positive margins")
    print(f"Minimum (sigma_k/m,L,margin): {minimum}")
    print(f"CSV: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    main(parser.parse_args().output)
