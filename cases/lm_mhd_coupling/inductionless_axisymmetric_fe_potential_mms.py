#!/usr/bin/env python3
"""Axisymmetric Q1 FE inductionless potential MMS.

This standalone reference assembles the axisymmetric weak form

    int R sigma grad(w).grad(phi) dR dZ = int R sigma grad(w).(u x B) dR dZ

with all-insulating natural current closure and a nodal mean-value gauge. It is
a metric/sign prototype only and does not run MUG or TokaMaker.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import smplotlib  # noqa: F401  (applies the style on import)
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from inductionless_fe_potential_mms import (  # noqa: E402
    observed_orders,
    q1_shapes,
    solve_with_nodal_mean_gauge,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--sizes", default="9,17,25,33,49", help="comma-separated grid sizes")
    return parser.parse_args()


def conductivity(r: float | np.ndarray, z: float | np.ndarray) -> float | np.ndarray:
    return 1.0 + 0.15 * (r - 1.0) + 0.10 * z


def phi_exact_raw(r: np.ndarray, z: np.ndarray) -> np.ndarray:
    return 0.015 * ((r - 1.5) ** 2 + (z - 0.5) ** 2)


def grad_phi_exact(r: float | np.ndarray, z: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    return 0.03 * (r - 1.5), 0.03 * (z - 0.5)


def target_current(r: float | np.ndarray, z: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    scale = 0.03
    psi_z = scale * (r - 1.0) * (2.0 - r) * (1.0 - 2.0 * z)
    psi_r = scale * (3.0 - 2.0 * r) * z * (1.0 - z)
    return psi_z / r, -psi_r / r


def motional_emf(r: float | np.ndarray, z: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    sigma = conductivity(r, z)
    grad_r, grad_z = grad_phi_exact(r, z)
    jr, jz = target_current(r, z)
    return grad_r + jr / sigma, grad_z + jz / sigma


def assemble_axisymmetric_system(r: np.ndarray, z: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    nr = r.size
    nz = z.size
    size = nr * nz
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    rhs = np.zeros(size, dtype=np.float64)
    gauss = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))

    def node(iz: int, ir: int) -> int:
        return iz * nr + ir

    for ez in range(nz - 1):
        for er in range(nr - 1):
            r0 = r[er]
            r1 = r[er + 1]
            z0 = z[ez]
            z1 = z[ez + 1]
            hr = r1 - r0
            hz = z1 - z0
            det_j = hr * hz / 4.0
            local_nodes = [node(ez, er), node(ez, er + 1), node(ez + 1, er + 1), node(ez + 1, er)]
            local_k = np.zeros((4, 4), dtype=np.float64)
            local_rhs = np.zeros(4, dtype=np.float64)
            for xi in gauss:
                for eta in gauss:
                    _, d_dxi, d_deta = q1_shapes(xi, eta)
                    rq = 0.5 * (r0 + r1) + 0.5 * hr * xi
                    zq = 0.5 * (z0 + z1) + 0.5 * hz * eta
                    sigma = float(conductivity(rq, zq))
                    emf_r, emf_z = motional_emf(rq, zq)
                    d_dr = d_dxi * (2.0 / hr)
                    d_dz = d_deta * (2.0 / hz)
                    weight = rq * sigma * det_j
                    for i in range(4):
                        local_rhs[i] += weight * (d_dr[i] * emf_r + d_dz[i] * emf_z)
                        for j in range(4):
                            local_k[i, j] += weight * (d_dr[i] * d_dr[j] + d_dz[i] * d_dz[j])
            for i, row in enumerate(local_nodes):
                rhs[row] += local_rhs[i]
                for j, col in enumerate(local_nodes):
                    rows.append(row)
                    cols.append(col)
                    data.append(local_k[i, j])
    return sp.csr_matrix((data, (rows, cols)), shape=(size, size)), rhs


def axisymmetric_divergence(current_r: np.ndarray, current_z: np.ndarray, r: np.ndarray, z: np.ndarray) -> np.ndarray:
    r_grid, _ = np.meshgrid(r, z, indexing="xy")
    d_rjr_dr = np.gradient(r_grid * current_r, r, axis=1, edge_order=2)
    d_jz_dz = np.gradient(current_z, z, axis=0, edge_order=2)
    return d_rjr_dr / r_grid + d_jz_dz


def solve_case(n: int) -> dict[str, np.ndarray | float]:
    r = np.linspace(1.0, 2.0, n)
    z = np.linspace(0.0, 1.0, n)
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    matrix, rhs = assemble_axisymmetric_system(r, z)
    phi = solve_with_nodal_mean_gauge(matrix, rhs).reshape(n, n)
    phi_exact = phi_exact_raw(r_grid, z_grid)
    phi_exact -= np.mean(phi_exact)
    error = phi - phi_exact
    sigma = conductivity(r_grid, z_grid)
    emf_r, emf_z = motional_emf(r_grid, z_grid)
    grad_r = np.gradient(phi, r, axis=1, edge_order=2)
    grad_z = np.gradient(phi, z, axis=0, edge_order=2)
    current_r = sigma * (-grad_r + emf_r)
    current_z = sigma * (-grad_z + emf_z)
    residual = axisymmetric_divergence(current_r, current_z, r, z)
    wall_current = max(
        float(np.max(np.abs(current_r[:, 0]))),
        float(np.max(np.abs(current_r[:, -1]))),
        float(np.max(np.abs(current_z[0, :]))),
        float(np.max(np.abs(current_z[-1, :]))),
    )
    return {
        "r": r,
        "z": z,
        "r_grid": r_grid,
        "z_grid": z_grid,
        "conductivity": sigma,
        "phi": phi,
        "phi_exact": phi_exact,
        "potential_error": error,
        "current_r": current_r,
        "current_z": current_z,
        "residual": residual,
        "h": float(1.0 / (n - 1)),
        "max_error": float(np.max(np.abs(error))),
        "rms_error": float(np.sqrt(np.mean(error**2))),
        "max_axisym_divJ": float(np.max(np.abs(residual[2:-2, 2:-2]))),
        "max_wall_current": wall_current,
    }


def write_metrics(outdir: Path, cases: list[dict[str, np.ndarray | float]]) -> None:
    errors = [float(case["max_error"]) for case in cases]
    hs = [float(case["h"]) for case in cases]
    orders = observed_orders(errors, hs)
    with (outdir / "inductionless_axisymmetric_fe_potential_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["n", "h", "max_error", "rms_error", "observed_order", "max_axisym_divJ", "max_wall_current"])
        for case, order in zip(cases, orders):
            writer.writerow(
                [
                    int(case["r"].size),
                    f"{float(case['h']):.12e}",
                    f"{float(case['max_error']):.12e}",
                    f"{float(case['rms_error']):.12e}",
                    "" if order is None else f"{order:.6f}",
                    f"{float(case['max_axisym_divJ']):.12e}",
                    f"{float(case['max_wall_current']):.12e}",
                ]
            )


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    r_grid = case["r_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    panels = [
        (case["conductivity"], r"$\sigma(R,Z)$", "viridis"),
        (case["phi_exact"], r"exact $\phi$", "magma"),
        (case["phi"], r"axisymmetric FE $\phi$", "magma"),
        (case["potential_error"], r"$\phi$ error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(r_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("R")
        ax.set_ylabel("Z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Axisymmetric Q1 weak-form potential MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_axisymmetric_fe_potential_mms_solution.png", dpi=180)
    plt.close(fig)


def plot_convergence(outdir: Path, cases: list[dict[str, np.ndarray | float]]) -> None:
    hs = np.array([float(case["h"]) for case in cases])
    errors = np.array([float(case["max_error"]) for case in cases])
    reference = errors[-1] * (hs / hs[-1]) ** 2
    fig, ax = plt.subplots(figsize=(6.3, 4.4))
    ax.loglog(hs, errors, "o-", label="axisymmetric Q1 max error")
    ax.loglog(hs, reference, "--", label=r"$O(h^2)$ reference")
    ax.invert_xaxis()
    ax.set_xlabel("h")
    ax.set_ylabel(r"max $|\phi-\phi_{exact}|$")
    ax.set_title("Axisymmetric weak-form convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_axisymmetric_fe_potential_mms_convergence.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, cases: list[dict[str, np.ndarray | float]]) -> None:
    errors = [float(case["max_error"]) for case in cases]
    hs = [float(case["h"]) for case in cases]
    orders = observed_orders(errors, hs)
    finest = cases[-1]
    lines = [
        "# Axisymmetric Q1 FE inductionless potential MMS",
        "",
        "This standalone Python reference assembles the weak form with the",
        "axisymmetric R metric. It does not run MUG or TokaMaker.",
        "",
        "| n | max error | observed order |",
        "|---:|---:|---:|",
    ]
    for case, order in zip(cases, orders):
        order_text = "" if order is None else f"`{order:.3f}`"
        lines.append(f"| {int(case['r'].size)} | `{float(case['max_error']):.3e}` | {order_text} |")
    lines.extend(
        [
            "",
            "Finest-grid nodal postprocessing:",
            f"- max axisymmetric |div J| = {float(finest['max_axisym_divJ']):.3e}.",
            f"- max wall-normal current = {float(finest['max_wall_current']):.3e}.",
            "",
            "Red-team interpretation:",
            "- This checks only the scalar axisymmetric R-weighted weak form.",
            "- It does not include toroidal electric field, tensor conductivity,",
            "  TokaMaker feedback, or MUG field extraction.",
            "- Current diagnostics are structured-grid nodal postprocessing, not",
            "  conservative FE flux recovery.",
        ]
    )
    (outdir / "inductionless_axisymmetric_fe_potential_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sizes = [int(item) for item in args.sizes.split(",") if item.strip()]
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cases = [solve_case(size) for size in sizes]
    write_metrics(outdir, cases)
    plot_solution(outdir, cases[-1])
    plot_convergence(outdir, cases)
    write_summary(outdir, cases)
    print(f"Wrote axisymmetric Q1 FE potential MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
