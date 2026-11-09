#!/usr/bin/env python3
"""Command-line entry point for the published spin-identifiability computations.

The full workflow regenerates the primary alpha=0 numerical tables. The
high-resolution adversarial search is optional and does not certify separation.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
FROZEN = ROOT / "theory_numerics" / "results" / "alpha0"
REPRODUCED = ROOT / "theory_numerics" / "reproduced_alpha0"


def execute(name: str, *args: str) -> None:
    command = [sys.executable, str(ROOT / "theory_numerics" / name), *args]
    print("$", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def smoke() -> None:
    classification = pd.read_csv(FROZEN / "classification_99_alpha0.csv")
    counts = classification["status"].value_counts().to_dict()
    assert counts == {"certified": 63, "overlap": 33, "unresolved": 3}, counts
    assert int(classification["componentwise_certified"].sum()) == 46
    frozen = pd.read_csv(FROZEN / "alpha0_48segments.csv")
    assert len(frozen) == 4
    assert frozen["all_slope_bounds_ok"].all()
    assert frozen["W_C_n20001"].iloc[0] < 0
    assert (frozen["W_C_n20001"].iloc[1:] > 0).all()
    execute("verify_48segment_quadrature.py")
    print("PASS frozen alpha=0 results and independent 48-segment quadrature")


def full() -> None:
    for script in (
        "spin_alpha0_compute.py",
        "spin_alpha0_qcqp.py",
        "spin_alpha0_power.py",
        "spin_alpha0_figures.py",
        "verify_alpha0_release.py",
    ):
        execute(script)
    smoke()
    print("PASS primary alpha=0 numerical reproduction")


def highres() -> None:
    REPRODUCED.mkdir(parents=True, exist_ok=True)
    execute("alpha0_48segment_audit.py", "--out", str(REPRODUCED), "--starts", "200")
    new = pd.read_csv(REPRODUCED / "alpha0_48segments.csv")
    assert len(new) == 4 and new["all_slope_bounds_ok"].all()
    assert new["W_C_n20001"].iloc[0] < 0
    print("PASS optional 48-segment search (positive values do not prove separation)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full", "highres"), default="smoke")
    args = parser.parse_args()
    {"smoke": smoke, "full": full, "highres": highres}[args.mode]()
