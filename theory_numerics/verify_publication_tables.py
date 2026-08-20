#!/usr/bin/env python3
"""Verify the archived finite-bandwidth publication tables and status counts."""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"

def main() -> None:
    s1 = pd.read_csv(RESULTS / "S1_slope_bandwidth_certificates_monotone.csv")
    counts = s1["status_monotone"].value_counts().to_dict()
    expected = {"guaranteed": 46, "explicit_overlap": 32, "unresolved": 21}
    assert counts == expected, (counts, expected)
    assert len(s1) == 99
    s2 = pd.read_csv(RESULTS / "S2_covariance_shape_scan.csv")
    ref = s2[(s2.r_sigmaC_over_sigmaM == 1) & (s2.rho == 0)].iloc[0]
    assert abs(ref.rigorous_lower_bound - 0.3682) < 5e-5
    s3 = pd.read_csv(RESULTS / "S3_precision_based_recovery.csv")
    assert abs(float(s3.loc[s3.tau == 0.50, "spin1_power_at_5pct_false_positive"].iloc[0]) - 0.881) < 5e-4
    assert abs(float(s3.loc[s3.tau == 0.15, "roc_auc"].iloc[0]) - 0.9997) < 5e-5
    summary = json.loads((RESULTS / "finite_bandwidth_publication_summary.json").read_text())
    assert summary["finite_bandwidth_grid"]["points"] == 99
    print("finite_bandwidth_publication_tables_ok")
    print(counts)

if __name__ == "__main__":
    main()
