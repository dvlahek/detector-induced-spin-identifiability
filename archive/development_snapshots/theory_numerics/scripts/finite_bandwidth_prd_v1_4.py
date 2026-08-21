#!/usr/bin/env python3
"""Reproduce the PRD v1.4 finite-bandwidth protected-sector classification.

The calculation separates three logically distinct statements.

1. Whole-Lipschitz certificate. Positivity of the magnetic detector measure
   implies that an origin overlap requires g_M=0 and G_M(x)=0. The remaining
   charge response is bounded below pointwise over the full Lipschitz class.
2. Finite-basis global check. In an n-segment piecewise-linear basis the charge
   response is a box-constrained QCQP. It is linear in the G1 slopes, so the
   global minimum is found by enumerating their 2**n vertices and solving one
   convex box QP in the GQ slopes at each vertex.
3. Degree-2 SOS certificate. A directly verifiable quadratic identity is
   optimized for the unreduced bilinear QCQP. The saved multipliers and PSD
   residuals make every successful certificate independently checkable.

The script also constructs the 32-segment explicit overlap at
(sigma_k/m,L)=(0.15,8), audits the two remaining unresolved points, and writes
the version-neutral 64/33/2 status map used in the PRD package.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize


SIGMAS = (0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
SLOPES = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0)
SEED = 20260820


@dataclass(frozen=True)
class SpectralData:
    k: np.ndarray
    x: np.ndarray
    mu_c: np.ndarray


@dataclass(frozen=True)
class ChargeQCQP:
    c0: float
    a_u: np.ndarray
    a_v: np.ndarray
    h_uv: np.ndarray
    q_vv: np.ndarray
    xmax: float


def spectral_data(sigma: float, grid_n: int = 2601,
                  kmax_sigma: float = 9.0) -> SpectralData:
    """Normalized conserved-charge measure used in the manuscript."""
    kmax = max(kmax_sigma * sigma, 1.0e-4)
    k = np.linspace(-kmax, kmax, grid_n)
    e = np.sqrt(1.0 + k * k)
    x = 2.0 * (e - 1.0)
    eta = x / 4.0
    w0 = np.exp(-0.5 * (k / sigma) ** 2) / (2.0 * e)
    conserved = (k / (e + 1.0)) ** 2
    w_c = w0 * conserved * eta
    norm = np.trapezoid(w_c, k)
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("non-positive charge normalization")
    return SpectralData(k=k, x=x, mu_c=w_c / norm)


def integrate(sd: SpectralData, values: np.ndarray) -> np.ndarray:
    """Integrate arrays whose first axis is the spectral grid."""
    shape = (len(sd.k),) + (1,) * (values.ndim - 1)
    return np.trapezoid(sd.mu_c.reshape(shape) * values, sd.k, axis=0)


def segment_basis(x: np.ndarray, xmax: float, n_segments: int) -> np.ndarray:
    """Integrated-slope basis phi_i(x) for a uniform piecewise-linear grid."""
    knots = np.linspace(0.0, xmax, n_segments + 1)
    delta = knots[1] - knots[0]
    return np.column_stack([
        np.clip(x - knots[i], 0.0, delta) for i in range(n_segments)
    ])


def charge_qcqp(sigma: float, n_segments: int) -> ChargeQCQP:
    sd = spectral_data(sigma)
    xmax = float(np.max(sd.x))
    phi = segment_basis(sd.x, xmax, n_segments)
    outer = phi[:, :, None] * phi[:, None, :]
    c0 = float(integrate(sd, 1.0 + sd.x / 4.0))
    a_u = np.asarray(integrate(sd, phi), dtype=float)
    a_v = np.asarray(integrate(sd, phi * (1.0 + sd.x[:, None] / 2.0)),
                     dtype=float)
    h_uv = np.asarray(integrate(sd, outer), dtype=float)
    q_vv = np.asarray(integrate(sd, outer * (sd.x / 4.0)[:, None, None]),
                      dtype=float)
    q_vv = 0.5 * (q_vv + q_vv.T)
    return ChargeQCQP(c0=c0, a_u=a_u, a_v=a_v, h_uv=h_uv,
                      q_vv=q_vv, xmax=xmax)


def charge_value(q: ChargeQCQP, u: np.ndarray, v: np.ndarray) -> float:
    return float(q.c0 + q.a_u @ u + q.a_v @ v
                 + u @ q.h_uv @ v + v @ q.q_vv @ v)


def pointwise_charge_envelope(x: np.ndarray, slope_bound: float) -> np.ndarray:
    """Exact box minimum h_L(x) used for the whole-Lipschitz certificate.

    For a=G1(x) and b=GQ(x), the anchored Lipschitz condition gives
    a,b in [1-Lx,1+Lx]. For fixed b the minimum over a is attained on an
    endpoint. Minimizing the resulting convex quadratic over b therefore only
    requires the two a endpoints.
    """
    lo = 1.0 - slope_bound * x
    hi = 1.0 + slope_bound * x
    out = np.empty_like(x)
    zero = x <= 1.0e-15
    out[zero] = 1.0
    nz = ~zero
    xn = x[nz]
    lon = lo[nz]
    hin = hi[nz]
    candidates = []
    for a in (lon, hin):
        b_star = np.clip(-2.0 * a / xn, lon, hin)
        candidates.append(a * b_star + 0.25 * xn * b_star * b_star)
    out[nz] = np.minimum(candidates[0], candidates[1])
    return out


def whole_lipschitz_lower_bound(sigma: float, slope_bound: float) -> float:
    sd = spectral_data(sigma, grid_n=20001, kmax_sigma=12.0)
    h = pointwise_charge_envelope(sd.x, slope_bound)
    return float(integrate(sd, h))


def solve_convex_v(q: ChargeQCQP, u: np.ndarray, slope_bound: float,
                   start: np.ndarray | None = None) -> tuple[float, np.ndarray, float]:
    """Solve the convex v-subproblem and return value, solution, KKT residual."""
    linear = q.a_v + q.h_uv.T @ u

    def fun(v: np.ndarray) -> float:
        return float(q.c0 + q.a_u @ u + linear @ v + v @ q.q_vv @ v)

    def jac(v: np.ndarray) -> np.ndarray:
        return linear + 2.0 * q.q_vv @ v

    if start is None:
        start = np.zeros_like(u)
    result = minimize(fun, np.clip(start, -slope_bound, slope_bound), jac=jac,
                      method="L-BFGS-B",
                      bounds=[(-slope_bound, slope_bound)] * len(u),
                      options={"ftol": 1.0e-15, "gtol": 1.0e-12,
                               "maxiter": 3000, "maxls": 80})
    v = np.asarray(result.x, dtype=float)
    grad = jac(v)
    projected = v - np.clip(v - grad, -slope_bound, slope_bound)
    residual = float(np.max(np.abs(projected)))
    return fun(v), v, residual


def global_finite_basis_minimum(sigma: float, slope_bound: float,
                                n_segments: int = 8) -> dict:
    """Global QCQP minimum by 2**n vertex enumeration plus convex QPs."""
    q = charge_qcqp(sigma, n_segments)
    best = (np.inf, None, None, np.inf)
    for signs in itertools.product((-1.0, 1.0), repeat=n_segments):
        u = slope_bound * np.asarray(signs, dtype=float)
        value, v, residual = solve_convex_v(q, u, slope_bound)
        if value < best[0]:
            best = (value, u.copy(), v.copy(), residual)
    return {
        "minimum": float(best[0]),
        "u": best[1],
        "v": best[2],
        "projected_kkt_residual": float(best[3]),
        "n_vertices": 2 ** n_segments,
        "n_segments": n_segments,
    }


def alternating_overlap_profile(sigma: float = 0.15, slope_bound: float = 8.0,
                                n_segments: int = 32,
                                starts: int = 512,
                                require_negative: bool = True) -> dict:
    """Search for a negative charge profile by exact block minimization.

    If ``require_negative`` is false, the best profile is returned even when
    its value remains positive. Such a result documents an adversarial search
    but is not promoted to a whole-class certificate.
    """
    q = charge_qcqp(sigma, n_segments)
    rng = np.random.default_rng(SEED)
    best = (np.inf, None, None, np.inf)
    initial_v = [np.zeros(n_segments),
                 np.full(n_segments, slope_bound),
                 np.full(n_segments, -slope_bound)]
    initial_v.extend(rng.uniform(-slope_bound, slope_bound, size=(starts, n_segments)))

    for v0 in initial_v:
        v = np.asarray(v0, dtype=float).copy()
        u_prev = None
        residual = np.inf
        for _ in range(200):
            coeff_u = q.a_u + q.h_uv @ v
            u = -slope_bound * np.where(coeff_u >= 0.0, 1.0, -1.0)
            value, v_new, residual = solve_convex_v(q, u, slope_bound, start=v)
            if u_prev is not None and np.array_equal(u, u_prev) \
                    and np.max(np.abs(v_new - v)) < 1.0e-11:
                v = v_new
                break
            u_prev = u.copy()
            v = v_new
        value = charge_value(q, u, v)
        if value < best[0]:
            best = (value, u.copy(), v.copy(), residual)

    if not np.isfinite(best[0]):
        raise RuntimeError("finite overlap-search profile not found")
    if best[0] >= 0.0 and require_negative:
        raise RuntimeError(f"negative overlap profile not found; best={best[0]:.12g}")

    target_u = best[1]
    target_v = best[2]

    def scaled_value(t: float) -> float:
        return charge_value(q, t * target_u, t * target_v)

    payload = {
        "sigma_k_over_m": sigma,
        "slope_bound_L": slope_bound,
        "n_segments": n_segments,
        "random_starts": starts,
        "xmax": q.xmax,
        "target_value_t0": scaled_value(0.0),
        "target_value_t1": scaled_value(1.0),
        "negative_profile_found": bool(best[0] < 0.0),
        "target_u": target_u,
        "target_v": target_v,
        "projected_kkt_residual": float(best[3]),
    }
    if best[0] >= 0.0:
        return payload

    root = float(brentq(scaled_value, 0.0, 1.0, xtol=5.0e-15,
                        rtol=5.0e-15, maxiter=300))
    u_root = root * target_u
    v_root = root * target_v
    payload.update({
        "root_t": root,
        "root_value": scaled_value(root),
        "max_abs_slope_at_root": float(max(np.max(np.abs(u_root)),
                                               np.max(np.abs(v_root)))),
        "root_u": u_root,
        "root_v": v_root,
    })
    return payload


def full_degree2_sos(q: ChargeQCQP, slope_bound: float) -> dict:
    """Optimize and verify the degree-2 SOS/S-procedure lower bound.

    For z=(u,v), write q(z)=c+l^T z+z^T R z. The certificate is

      q(z)-gamma = [1,z]^T S [1,z]
                   + sum_i lambda_i (L^2-z_i^2),

    with S positive semidefinite and lambda_i nonnegative. Eliminating S and
    gamma reduces the SDP dual to a smooth convex minimization over lambda on
    the domain R+diag(lambda) positive definite.
    """
    n = len(q.a_u)
    linear = np.concatenate([q.a_u, q.a_v])
    rmat = np.block([
        [np.zeros((n, n)), 0.5 * q.h_uv],
        [0.5 * q.h_uv.T, q.q_vv],
    ])
    rmat = 0.5 * (rmat + rmat.T)
    dim = 2 * n
    shift = max(1.0e-6, -float(np.linalg.eigvalsh(rmat)[0]) + 1.0e-4)
    lam0 = np.full(dim, shift)

    def objective(lam: np.ndarray) -> tuple[float, np.ndarray]:
        amat = rmat + np.diag(lam)
        eigmin = float(np.linalg.eigvalsh(amat)[0])
        if eigmin <= 1.0e-13:
            penalty = 1.0e12 + 1.0e12 * (1.0e-13 - eigmin)
            return penalty, np.full(dim, -1.0e6)
        y = np.linalg.solve(amat, linear)
        value = slope_bound ** 2 * np.sum(lam) + 0.25 * linear @ y - q.c0
        grad = slope_bound ** 2 - 0.25 * y * y
        return float(value), np.asarray(grad, dtype=float)

    result = minimize(lambda x: objective(x)[0], lam0,
                      jac=lambda x: objective(x)[1], method="L-BFGS-B",
                      bounds=[(0.0, None)] * dim,
                      options={"ftol": 1.0e-14, "gtol": 1.0e-10,
                               "maxiter": 10000, "maxls": 100})
    lam = np.maximum(np.asarray(result.x, dtype=float), 0.0)
    amat = rmat + np.diag(lam)
    y = np.linalg.solve(amat, linear)
    gamma_opt = float(q.c0 - slope_bound ** 2 * np.sum(lam)
                      - 0.25 * linear @ y)

    # A small downward rounding makes the saved matrix a numerical certificate.
    gamma = gamma_opt - 1.0e-10 * max(1.0, abs(gamma_opt))
    smat = np.empty((dim + 1, dim + 1), dtype=float)
    smat[0, 0] = q.c0 - gamma - slope_bound ** 2 * np.sum(lam)
    smat[0, 1:] = 0.5 * linear
    smat[1:, 0] = 0.5 * linear
    smat[1:, 1:] = amat
    eigs = np.linalg.eigvalsh(0.5 * (smat + smat.T))

    # Coefficient reconstruction is exact by construction; report floating
    # residuals explicitly so certificates can be audited without trusting the
    # optimizer status flag.
    recovered_c = smat[0, 0] + slope_bound ** 2 * np.sum(lam) + gamma
    recovered_l = 2.0 * smat[0, 1:]
    recovered_r = smat[1:, 1:] - np.diag(lam)
    coefficient_residual = max(
        abs(recovered_c - q.c0),
        float(np.max(np.abs(recovered_l - linear))),
        float(np.max(np.abs(recovered_r - rmat))),
    )
    return {
        "gamma_optimized": gamma_opt,
        "gamma_verified": gamma,
        "lambda": lam,
        "minimum_psd_eigenvalue": float(eigs[0]),
        "coefficient_residual": float(coefficient_residual),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
    }


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def plot_status_map(df: pd.DataFrame, output: Path) -> None:
    code = {"guaranteed": 1, "unresolved": 0, "explicit_overlap": -1}
    pivot = df.assign(code=df.status_prd_v1_4.map(code)).pivot(
        index="L", columns="sigma_k_over_m", values="code"
    ).reindex(index=SLOPES, columns=SIGMAS)
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cmap = ListedColormap(["#d95f5f", "#d9d9d9", "#58a66a"])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
    image = ax.imshow(pivot.values, origin="lower", aspect="auto",
                      cmap=cmap, norm=norm)
    ax.set_xticks(range(len(SIGMAS)), [f"{s:g}" for s in SIGMAS])
    ax.set_yticks(range(len(SLOPES)), [f"{ell:g}" for ell in SLOPES])
    ax.set_xlabel(r"kernel width $\sigma_k/m$")
    ax.set_ylabel(r"Lipschitz bound $L$")
    ax.set_title("Finite-bandwidth protected-sector classification")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", markersize=10,
                          markerfacecolor=color, markeredgecolor="none", label=label)
               for color, label in [("#58a66a", "whole-class certified"),
                                    ("#d9d9d9", "unresolved"),
                                    ("#d95f5f", "explicit overlap")]]
    ax.legend(handles=handles, loc="upper left", frameon=True)
    fig.colorbar(image, ax=ax, ticks=[]).remove()
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-status", type=Path, required=True,
                        help="v1.1 S1_slope_bandwidth_certificates_monotone.csv")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--overlap-starts", type=int, default=512)
    parser.add_argument("--unresolved-overlap-starts", type=int, default=1024)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    old = pd.read_csv(args.old_status)
    old = old.rename(columns={"status_monotone": "status_v1_1",
                              "slope_bound_L": "L"})
    if "L" not in old:
        raise ValueError("old status table lacks L column")

    envelope_rows = []
    for sigma in SIGMAS:
        for ell in SLOPES:
            envelope_rows.append({
                "sigma_k_over_m": sigma,
                "L": ell,
                "whole_lipschitz_charge_lower_bound":
                    whole_lipschitz_lower_bound(sigma, ell),
            })
    envelope = pd.DataFrame(envelope_rows)
    merged = old.merge(envelope, on=["sigma_k_over_m", "L"], how="left",
                       validate="one_to_one")

    overlap = alternating_overlap_profile(starts=args.overlap_starts)
    with (args.outdir / "explicit_overlap_0p15_L8.json").open("w") as handle:
        json.dump(jsonable(overlap), handle, indent=2)

    def classify(row: pd.Series) -> str:
        if row.status_v1_1 == "guaranteed":
            return "guaranteed"
        if row.status_v1_1 == "explicit_overlap":
            return "explicit_overlap"
        if row.whole_lipschitz_charge_lower_bound > 1.0e-10:
            return "guaranteed"
        if abs(row.sigma_k_over_m - 0.15) < 1.0e-12 \
                and abs(row.L - 8.0) < 1.0e-12:
            return "explicit_overlap"
        return "unresolved"

    merged["status_prd_v1_4"] = merged.apply(classify, axis=1)
    merged.to_csv(args.outdir / "finite_bandwidth_status_prd_v1_4.csv", index=False)
    plot_status_map(merged, args.outdir / "finite_bandwidth_status_prd_v1_4.png")

    unresolved_searches = [
        alternating_overlap_profile(
            sigma=sigma,
            slope_bound=ell,
            n_segments=32,
            starts=args.unresolved_overlap_starts,
            require_negative=False,
        )
        for sigma, ell in ((0.08, 24.0), (0.10, 16.0))
    ]
    with (args.outdir / "unresolved_32segment_search_prd_v1_4.json").open("w") as handle:
        json.dump(jsonable(unresolved_searches), handle, indent=2)

    old_unresolved = merged[merged.status_v1_1 == "unresolved"].copy()
    qcqp_rows = []
    sos_payload = []
    for row in old_unresolved.itertuples(index=False):
        qcqp = global_finite_basis_minimum(row.sigma_k_over_m, row.L, 8)
        sos = full_degree2_sos(charge_qcqp(row.sigma_k_over_m, 8), row.L)
        qcqp_rows.append({
            "sigma_k_over_m": row.sigma_k_over_m,
            "L": row.L,
            "global_8segment_minimum": qcqp["minimum"],
            "projected_kkt_residual": qcqp["projected_kkt_residual"],
            "vertices_checked": qcqp["n_vertices"],
            "whole_lipschitz_charge_lower_bound":
                row.whole_lipschitz_charge_lower_bound,
            "status_prd_v1_4": row.status_prd_v1_4,
            "sos_gamma_verified": sos["gamma_verified"],
            "sos_minimum_psd_eigenvalue": sos["minimum_psd_eigenvalue"],
            "sos_coefficient_residual": sos["coefficient_residual"],
            "sos_positive_certificate": bool(
                sos["gamma_verified"] > 0.0
                and sos["minimum_psd_eigenvalue"] >= -1.0e-9
                and sos["coefficient_residual"] <= 1.0e-10
            ),
        })
        sos_payload.append({
            "sigma_k_over_m": row.sigma_k_over_m,
            "L": row.L,
            "qcqp": qcqp,
            "sos": sos,
        })

    qcqp_df = pd.DataFrame(qcqp_rows)
    qcqp_df.to_csv(args.outdir / "qcqp_sos_unresolved_prd_v1_4.csv", index=False)
    with (args.outdir / "qcqp_sos_certificates_prd_v1_4.json").open("w") as handle:
        json.dump(jsonable(sos_payload), handle, indent=2)

    counts = merged.status_prd_v1_4.value_counts().to_dict()
    unresolved = merged[merged.status_prd_v1_4 == "unresolved"][
        ["sigma_k_over_m", "L"]
    ].to_dict(orient="records")
    summary = {
        "status_counts": counts,
        "unresolved_points": unresolved,
        "old_unresolved_points": int(len(old_unresolved)),
        "old_unresolved_whole_class_certified": int(
            np.sum((old_unresolved.whole_lipschitz_charge_lower_bound > 1.0e-10))
        ),
        "old_unresolved_sos_certified": int(qcqp_df.sos_positive_certificate.sum()),
        "all_old_unresolved_8segment_minima_positive": bool(
            np.all(qcqp_df.global_8segment_minimum > 0.0)
        ),
        "maximum_8segment_kkt_residual": float(
            qcqp_df.projected_kkt_residual.max()
        ),
        "explicit_overlap": {
            key: overlap[key] for key in (
                "sigma_k_over_m", "slope_bound_L", "n_segments",
                "target_value_t0", "target_value_t1", "root_t",
                "root_value", "max_abs_slope_at_root"
            )
        },
        "unresolved_32segment_search": [
            {
                key: search[key] for key in (
                    "sigma_k_over_m", "slope_bound_L", "n_segments",
                    "random_starts", "target_value_t1",
                    "negative_profile_found", "projected_kkt_residual",
                )
            }
            for search in unresolved_searches
        ],
    }
    with (args.outdir / "finite_bandwidth_prd_v1_4_summary.json").open("w") as handle:
        json.dump(jsonable(summary), handle, indent=2)
    print(json.dumps(jsonable(summary), indent=2))


if __name__ == "__main__":
    main()
