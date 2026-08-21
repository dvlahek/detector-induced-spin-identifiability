#!/usr/bin/env python3
"""Fail-fast consistency checks for the publication v1.0 publication tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "theory_numerics" / "results" / "publication_v1_0"


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if not np.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError(
            f"{label}: got {actual:.16g}, expected {expected:.16g} "
            f"within {tolerance:.3g}"
        )


def main(data: Path) -> None:
    DATA = data
    status = pd.read_csv(DATA / "finite_bandwidth_status.csv")
    counts = status.status_publication.value_counts().to_dict()
    if counts != {"guaranteed": 64, "explicit_overlap": 33, "unresolved": 2}:
        raise AssertionError(f"unexpected phase-map counts: {counts}")

    diagnostics = pd.read_csv(DATA / "detector_measure_diagnostics.csv")
    row = diagnostics[np.isclose(diagnostics.sigma_k_over_m, 0.15)].iloc[0]
    close(float(row.mean_x_C), 0.10344928430838263, 2.0e-12,
          "charge-measure mean")
    close(float(row.normalization_M), 1.0, 5.0e-13,
          "magnetic normalization")
    close(float(row.normalization_C), 1.0, 5.0e-13,
          "charge normalization")

    coordinate = pd.read_csv(DATA / "coordinate_robustness.csv")
    coordinate_counts = coordinate.groupby("coordinate").origin_envelope_positive.sum()
    if coordinate_counts.to_dict() != {"G1": 64, "GC": 63}:
        raise AssertionError(f"unexpected coordinate counts: {coordinate_counts.to_dict()}")

    covariance = pd.read_csv(DATA / "whole_class_distance_covariance.csv")
    reference = covariance[
        np.isclose(covariance.sigmaC_over_sigmaM, 1.0)
        & np.isclose(covariance.rho, 0.0)
    ].iloc[0]
    close(float(reference.whole_class_whitened_lower_bound),
          0.40425535950546276, 2.0e-10, "reference whole-class distance")
    broad_charge = covariance[
        np.isclose(covariance.sigmaC_over_sigmaM, 4.0)
        & np.isclose(covariance.rho, 0.0)
    ].iloc[0]
    close(float(broad_charge.whole_class_whitened_lower_bound),
          0.15092614483672623, 2.0e-10, "broad-charge distance")

    convergence = pd.read_csv(DATA / "distance_quadrature_convergence.csv")
    for convention in ("G1", "GC"):
        values = convergence[
            convergence.coordinate == convention
        ].whole_class_whitened_lower_bound.to_numpy()
        if float(np.ptp(values)) > 1.0e-9:
            raise AssertionError(f"quadrature spread too large for {convention}")

    overlap = json.loads((DATA / "explicit_overlap_0p15_L8.json").read_text())
    close(abs(float(overlap["root_value"])), 0.0, 2.0e-14,
          "explicit overlap root")
    if float(overlap["max_abs_slope_at_root"]) >= 8.0:
        raise AssertionError("explicit overlap violates slope bound")

    bracket = json.loads(
        (DATA / "reference_distance_bracket.json").read_text())
    close(float(bracket["lower_bound"]), 0.40425535950546276,
          2.0e-10, "distance-bracket lower endpoint")
    close(float(bracket["explicit_profile_upper_bound"]),
          0.40426571752366414, 2.0e-10,
          "distance-bracket upper endpoint")
    if float(bracket["explicit_profile_upper_bound"]) < float(bracket["lower_bound"]):
        raise AssertionError("distance bracket is reversed")

    distance_grid = pd.read_csv(DATA / "whole_class_distance_grid.csv")
    if len(distance_grid) != 99:
        raise AssertionError(f"unexpected distance-grid size: {len(distance_grid)}")
    if int((distance_grid.whole_class_distance_lower_bound > 1.0e-10).sum()) != 64:
        raise AssertionError("distance grid does not reproduce 64 positive bounds")

    sachs = json.loads(
        (DATA / "sachs_overlap_witness_audit.json").read_text())
    close(float(sachs["max_abs_dGC_dx"]), 7.927784804388963,
          2.0e-10, "Sachs witness slope")
    if not bool(sachs["sachs_lipschitz_admissible"]):
        raise AssertionError("explicit overlap is not Sachs-admissible")

    print("Publication 1.0 checks passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    main(args.data)
