#!/usr/bin/env python3
"""Static nuisance-domain sensitivity for the protected spin certificate.

For the elementary spin-1 static fingerprint
    W_M(g) = g^2,  W_C(g) = 1-g,
compute the uniform single-channel and combined margins over representative
compact nuisance intervals K. The script also produces the centered-family
scan K_R=[2-R,2+R] used in the manuscript.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures" / "nuisance_domain_sensitivity.png"
CSV = ROOT / "results" / "nuisance_domain_sensitivity.csv"

# Unique real root of 2 g^3 + g - 1 = 0.
roots = np.roots([2.0, 0.0, 1.0, -1.0])
g_star = float(np.real(roots[np.isreal(roots)][0]))


def margins(a: float, b: float):
    if a > b:
        raise ValueError("Require a <= b")

    # Magnetic-only margin inf_K |W_M| = inf_K g^2.
    delta_m = 0.0 if a <= 0.0 <= b else min(a * a, b * b)

    # Charge-only margin inf_K |W_C| = inf_K |1-g|.
    delta_c = 0.0 if a <= 1.0 <= b else min(abs(1.0 - a), abs(1.0 - b))

    # Combined Euclidean margin. f(g)=g^4+(1-g)^2 is strictly convex.
    g_proj = min(max(g_star, a), b)
    delta_mc = float(np.sqrt(g_proj**4 + (1.0 - g_proj) ** 2))

    tau = 1 if (delta_m > 0.0 or delta_c > 0.0) else 2
    return delta_m, delta_c, delta_mc, g_proj, tau


intervals = [
    ("{2}", 2.0, 2.0),
    ("[1.9,2.1]", 1.9, 2.1),
    ("[1.75,2.25]", 1.75, 2.25),
    ("[1.5,2.5]", 1.5, 2.5),
    ("[1,3]", 1.0, 3.0),
    ("[0.5,3.5]", 0.5, 3.5),
    ("[0,4]", 0.0, 4.0),
    ("[0,1]", 0.0, 1.0),
    ("[-4,4]", -4.0, 4.0),
]

rows = []
for label, a, b in intervals:
    dm, dc, dmc, gp, tau = margins(a, b)
    rows.append(
        {
            "K": label,
            "g_min": a,
            "g_max": b,
            "magnetic_margin": dm,
            "charge_margin": dc,
            "combined_margin": dmc,
            "combined_minimizer_g": gp,
            "tau_protected": tau,
        }
    )

CSV.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(CSV, index=False)

R = np.linspace(0.0, 2.5, 501)
dm_curve, dc_curve, dmc_curve = [], [], []
for r in R:
    dm, dc, dmc, _, _ = margins(2.0 - r, 2.0 + r)
    dm_curve.append(dm)
    dc_curve.append(dc)
    dmc_curve.append(dmc)

FIG.parent.mkdir(parents=True, exist_ok=True)
fig, ax = plt.subplots(figsize=(7.6, 4.8))
ax.plot(R, dm_curve, label=r"magnetic-only $\Delta_M$")
ax.plot(R, dc_curve, label=r"charge-only $\Delta_C$")
ax.plot(R, dmc_curve, label=r"combined $\Delta_{MC}$")
ax.axvline(1.0, linestyle="--", linewidth=1.0, label=r"$K_R$ reaches $g=1$")
ax.axvline(2.0, linestyle=":", linewidth=1.2, label=r"$K_R$ reaches $g=0$")
ax.text(0.18, 0.13, r"$\tau(K_R)=1$", transform=ax.transAxes)
ax.text(0.83, 0.13, r"$\tau(K_R)=2$", transform=ax.transAxes)
ax.set_xlabel(r"half-width $R$ for $K_R=[2-R,2+R]$")
ax.set_ylabel("uniform protected-response margin")
ax.set_xlim(0.0, 2.5)
ax.set_ylim(bottom=-0.03)
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(FIG, dpi=240)
plt.close(fig)

# Exact publication checkpoints.
assert abs(g_star - 0.5897545123) < 1e-10
assert abs(np.sqrt(g_star**4 + (1.0 - g_star) ** 2) - 0.5378414487) < 1e-10
assert margins(1.5, 2.5)[4] == 1
assert margins(0.5, 3.5)[4] == 1
assert margins(0.0, 4.0)[4] == 2
assert margins(-4.0, 4.0)[4] == 2

print(f"g_star={g_star:.10f}")
print(f"combined_global_min={np.sqrt(g_star**4 + (1.0 - g_star) ** 2):.10f}")
print(f"csv={CSV}")
print(f"figure={FIG}")
print("nuisance_domain_sensitivity_ok")
