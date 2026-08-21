#!/usr/bin/env python3
"""Apply monotone closure to a finite-bandwidth status table.

The Lipschitz classes are nested in L: if an explicit overlap is found for a
slope bound L0 at fixed sigma_k/m, the same counterexample is admissible for
all L >= L0. Guaranteed points are left unchanged because they come from an
independent analytic lower-bound certificate. Remaining points are unresolved.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def apply_monotone_closure(df: pd.DataFrame) -> pd.DataFrame:
    required = {"sigma_k_over_m", "L", "status"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    out = df.copy()
    for sigma, idx in out.groupby("sigma_k_over_m").groups.items():
        sub = out.loc[idx].sort_values("L")
        overlaps = sub.loc[sub["status"] == "explicit_overlap", "L"]
        if len(overlaps):
            first = overlaps.min()
            mask = (out["sigma_k_over_m"] == sigma) & (out["L"] >= first) & (out["status"] != "guaranteed")
            out.loc[mask, "status"] = "explicit_overlap"
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    df = pd.read_csv(args.input)
    out = apply_monotone_closure(df)
    out.to_csv(args.output, index=False)
    print(out["status"].value_counts().to_string())


if __name__ == "__main__":
    main()
