#!/usr/bin/env python3
"""Certify the 33 alpha=0 protected-charge endpoint signs with interval arithmetic.

Run from the repository root:
    python theory_numerics/verify_overlap_sign_intervals.py

The verdict compares interval upper endpoints directly with zero. The
CSV/JSON decimal values are diagnostics, not inputs to the sign test.
"""
import argparse
import csv
import json
from decimal import Decimal
from pathlib import Path
from mpmath import iv

DEFAULT_CSV = Path("theory_numerics/results/alpha0/explicit_overlaps_alpha0.csv")
DEFAULT_OUT = Path("theory_numerics/reproduced_alpha0/interval_endpoint_signs_alpha0.json")


def interval_min(a, b):
    return iv.mpf([a.a if a.a < b.a else b.a, a.b if a.b < b.b else b.b])


def certify(sigma_text, L_text, n=1000, cutoff=14):
    sigma_d = Decimal(str(sigma_text))
    L_d = Decimal(str(L_text))
    radius_d = sigma_d * Decimal(cutoff)
    step_d = radius_d / n
    sigma, L = iv.mpf(str(sigma_d)), iv.mpf(str(L_d))
    X = 2 * (iv.sqrt(1 + (9 * sigma) ** 2) - 1)
    H = X / 32
    integral = iv.mpf(0)
    for i in range(n):
        left = step_d * i
        right = left + step_d
        r = iv.mpf([str(left), str(right)])
        x = 2 * (iv.sqrt(1 + r*r) - 1)
        gq = 1 - L * interval_min(x, X)
        g1 = 1 + L * (interval_min(x, X) - 2 * interval_min(x, H))
        f = g1*gq + x*gq*gq/4
        w = iv.exp(-r*r/(2*sigma*sigma))*r*r*x/(2*iv.sqrt(1+r*r))
        integral += w*f*iv.mpf(str(step_d))

    # For r>=R: x<=2r, |G1|,|GQ|<=1+2Lr, and
    # w(r)<=r^3 exp[-r^2/(2sigma^2)].
    # Bound Gaussian moments J_p via integration by parts.
    R = iv.mpf(str(radius_d))
    exp_tail = iv.exp(-R*R/(2*sigma*sigma))
    j0_hi = sigma*sigma/R*exp_tail
    j1 = sigma*sigma*exp_tail
    j2_hi = sigma*sigma*R*exp_tail + sigma*sigma*j0_hi
    j3 = sigma*sigma*R**2*exp_tail + 2*sigma*sigma*j1
    j4_hi = sigma*sigma*R**3*exp_tail + 3*sigma*sigma*j2_hi
    j5 = sigma*sigma*R**4*exp_tail + 4*sigma*sigma*j3
    j6_hi = sigma*sigma*R**5*exp_tail + 5*sigma*sigma*j4_hi
    tail = (j3 + (4*L+iv.mpf('0.5'))*j4_hi
            + (4*L*L+2*L)*j5 + 2*L*L*j6_hi)
    bound = integral + iv.mpf([-tail.b, tail.b])
    return {
        "sigma_k_over_m": str(sigma_d), "L": str(L_d),
        "lower_interval": str(bound.a), "upper_interval": str(bound.b),
        "upper_report": float(bound.b),
        "tail_upper_report": float(tail.b),
        "negative_sign_proved_by_interval_comparison": bool(bound.b < iv.mpf(0)),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--subintervals", type=int, default=1000)
    p.add_argument("--cutoff", type=int, default=14)
    args = p.parse_args()
    if args.subintervals < 2 or args.cutoff < 9:
        p.error("Require >=2 subintervals and >=9 sigma cutoff")
    iv.dps = 35
    with args.csv.open(encoding="utf-8", newline="") as file:
        cases = list(csv.DictReader(file))
    if len(cases) != 33:
        raise RuntimeError(f"Expected 33 published cases, found {len(cases)}")
    rows = [certify(c["sigma_k_over_m"], c["L"],
                    args.subintervals, args.cutoff) for c in cases]
    if not all(r["negative_sign_proved_by_interval_comparison"] for r in rows):
        raise RuntimeError("At least one endpoint sign is not certified")
    worst = max(rows, key=lambda r: r["upper_report"])
    report = {
        "count": len(rows), "all_upper_intervals_strictly_negative": True,
        "interval_precision_digits": 35, "subintervals": args.subintervals,
        "cutoff_sigma": args.cutoff,
        "largest_upper_enclosure": worst, "results": rows,
        "scope": "1D alpha=0 protected-response endpoints, not full QFT outcome laws",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("INTERVAL_ENDPOINT_SIGN_33_PASS", worst["sigma_k_over_m"],
          worst["L"], worst["upper_interval"])
    print("Report:", args.out)


if __name__ == "__main__":
    main()
