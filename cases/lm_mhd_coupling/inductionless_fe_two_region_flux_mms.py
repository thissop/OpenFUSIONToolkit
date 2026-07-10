#!/usr/bin/env python3
"""Q1 FE two-region inductionless interface-flux MMS.

This standalone reference uses a conductivity jump aligned with element faces.
It solves the insulating weak form and recovers normal current on both sides of
the material interface. It is not a MUG implementation.
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
    fe_gradient_at,
    q1_shapes,
    solve_with_nodal_mean_gauge,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--sizes", default="17,33,49", help="comma-separated odd grid sizes")
    parser.add_argument("--sigma-left", type=float, default=1.0, help="left-region conductivity")
    parser.add_argument("--sigma-right", type=float, default=4.0, help="right-region conductivity")
    parser.add_argument("--flux", type=float, default=0.7, help="constant sigma grad(phi) flux")
    parser.add_argument("--current-scale", type=float, default=0.04, help="polynomial streamfunction amplitude")
    return parser.parse_args()


def sigma_region(x: float | np.ndarray, sigma_left: float, sigma_right: float) -> float | np.ndarray:
    return np.where(np.asarray(x) < 0.5, sigma_left, sigma_right)


def phi_exact_raw(
    x: np.ndarray,
    sigma_left: float,
    sigma_right: float,
    flux: float,
) -> np.ndarray:
    return np.where(
        x <= 0.5,
        flux * x / sigma_left,
        flux * 0.5 / sigma_left + flux * (x - 0.5) / sigma_right,
    )


def target_current(
    x: float | np.ndarray,
    z: float | np.ndarray,
    scale: float,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    return (
        scale * x * (1.0 - x) * (1.0 - 2.0 * z),
        -scale * (1.0 - 2.0 * x) * z * (1.0 - z),
    )


def motional_emf(
    x: float | np.ndarray,
    z: float | np.ndarray,
    sigma: float | np.ndarray,
    flux: float,
    current_scale: float,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    jx, jz = target_current(x, z, current_scale)
    return flux / sigma + jx / sigma, jz / sigma


def assemble_system(
    x: np.ndarray,
    z: np.ndarray,
    sigma_left: float,
    sigma_right: float,
    flux: float,
    current_scale: float,
) -> tuple[sp.csr_matrix, np.ndarray]:
    nx = x.size
    nz = z.size
    size = nx * nz
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    rhs = np.zeros(size, dtype=np.float64)
    gauss = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))

    def node(iz: int, ix: int) -> int:
        return iz * nx + ix

    for ez in range(nz - 1):
        for ex in range(nx - 1):
            x0 = x[ex]
            x1 = x[ex + 1]
            z0 = z[ez]
            z1 = z[ez + 1]
            hx = x1 - x0
            hz = z1 - z0
            det_j = hx * hz / 4.0
            local_nodes = [node(ez, ex), node(ez, ex + 1), node(ez + 1, ex + 1), node(ez + 1, ex)]
            local_k = np.zeros((4, 4), dtype=np.float64)
            local_rhs = np.zeros(4, dtype=np.float64)
            for xi in gauss:
                for eta in gauss:
                    _, d_dxi, d_deta = q1_shapes(xi, eta)
                    xq = 0.5 * (x0 + x1) + 0.5 * hx * xi
                    zq = 0.5 * (z0 + z1) + 0.5 * hz * eta
                    sigma = float(sigma_region(xq, sigma_left, sigma_right))
                    emf_x, emf_z = motional_emf(xq, zq, sigma, flux, current_scale)
                    d_dx = d_dxi * (2.0 / hx)
                    d_dz = d_deta * (2.0 / hz)
                    for i in range(4):
                        local_rhs[i] += sigma * (d_dx[i] * emf_x + d_dz[i] * emf_z) * det_j
                        for j in range(4):
                            local_k[i, j] += sigma * (d_dx[i] * d_dx[j] + d_dz[i] * d_dz[j]) * det_j
            for i, row in enumerate(local_nodes):
                rhs[row] += local_rhs[i]
                for j, col in enumerate(local_nodes):
                    rows.append(row)
                    cols.append(col)
                    data.append(local_k[i, j])
    return sp.csr_matrix((data, (rows, cols)), shape=(size, size)), rhs


def interface_fluxes(
    phi: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    sigma_left: float,
    sigma_right: float,
    flux: float,
    current_scale: float,
) -> dict[str, np.ndarray]:
    interface_idx = np.where(np.isclose(x, 0.5))[0]
    if interface_idx.size != 1:
        raise ValueError("x grid must contain the material interface at x=0.5")
    ix = int(interface_idx[0])
    gauss = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))
    z_values = []
    left_jx = []
    right_jx = []
    target_jx = []
    for ez in range(z.size - 1):
        z0 = z[ez]
        z1 = z[ez + 1]
        for eta in gauss:
            zq = 0.5 * (z0 + z1) + 0.5 * (z1 - z0) * eta
            grad_left, _ = fe_gradient_at(phi, x, z, ez, ix - 1, 1.0, eta)
            grad_right, _ = fe_gradient_at(phi, x, z, ez, ix, -1.0, eta)
            emf_left, _ = motional_emf(0.5, zq, sigma_left, flux, current_scale)
            emf_right, _ = motional_emf(0.5, zq, sigma_right, flux, current_scale)
            target, _ = target_current(0.5, zq, current_scale)
            z_values.append(zq)
            left_jx.append(sigma_left * (-grad_left + emf_left))
            right_jx.append(sigma_right * (-grad_right + emf_right))
            target_jx.append(target)
    return {
        "z": np.asarray(z_values),
        "left_jx": np.asarray(left_jx),
        "right_jx": np.asarray(right_jx),
        "target_jx": np.asarray(target_jx),
    }


def solve_case(args: argparse.Namespace, n: int) -> dict[str, np.ndarray | float | dict[str, np.ndarray]]:
    if n % 2 == 0:
        raise ValueError("FE two-region sizes must be odd so x=0.5 is a node line")
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.0, 1.0, n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    matrix, rhs = assemble_system(x, z, args.sigma_left, args.sigma_right, args.flux, args.current_scale)
    phi = solve_with_nodal_mean_gauge(matrix, rhs).reshape(n, n)
    phi_1d = phi_exact_raw(x, args.sigma_left, args.sigma_right, args.flux)
    phi_exact = np.broadcast_to(phi_1d, x_grid.shape).copy()
    phi_exact -= np.mean(phi_exact)
    sigma_grid = sigma_region(x_grid, args.sigma_left, args.sigma_right)
    fluxes = interface_fluxes(phi, x, z, args.sigma_left, args.sigma_right, args.flux, args.current_scale)
    interface_target_error = np.maximum(
        np.abs(fluxes["left_jx"] - fluxes["target_jx"]),
        np.abs(fluxes["right_jx"] - fluxes["target_jx"]),
    )
    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "conductivity": sigma_grid,
        "phi": phi,
        "phi_exact": phi_exact,
        "potential_error": phi - phi_exact,
        "interface_fluxes": fluxes,
        "h": float(1.0 / (n - 1)),
        "max_error": float(np.max(np.abs(phi - phi_exact))),
        "rms_error": float(np.sqrt(np.mean((phi - phi_exact) ** 2))),
        "max_interface_target_error": float(np.max(interface_target_error)),
        "max_interface_jump": float(np.max(np.abs(fluxes["left_jx"] - fluxes["right_jx"]))),
    }


def write_metrics(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, np.ndarray]]]) -> None:
    with (outdir / "inductionless_fe_two_region_flux_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "n",
                "h",
                "max_error",
                "rms_error",
                "max_interface_target_error",
                "max_interface_jump",
            ]
        )
        for case in cases:
            writer.writerow(
                [
                    int(case["x"].size),
                    f"{float(case['h']):.12e}",
                    f"{float(case['max_error']):.12e}",
                    f"{float(case['rms_error']):.12e}",
                    f"{float(case['max_interface_target_error']):.12e}",
                    f"{float(case['max_interface_jump']):.12e}",
                ]
            )


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float | dict[str, np.ndarray]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    panels = [
        (case["conductivity"], r"$\sigma$", "viridis"),
        (case["phi_exact"], r"exact $\phi$", "magma"),
        (case["phi"], r"FE $\phi$", "magma"),
        (case["potential_error"], r"$\phi$ error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.axvline(0.5, color="white", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Q1 FE two-region interface-flux MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_fe_two_region_flux_mms_solution.png", dpi=180)
    plt.close(fig)


def plot_interface(outdir: Path, case: dict[str, np.ndarray | float | dict[str, np.ndarray]]) -> None:
    fluxes = case["interface_fluxes"]
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.1), sharex=True)
    axes[0].plot(fluxes["z"], fluxes["target_jx"], label="target")
    axes[0].plot(fluxes["z"], fluxes["left_jx"], "--", label="left FE face")
    axes[0].plot(fluxes["z"], fluxes["right_jx"], ":", label="right FE face")
    axes[0].set_xlabel("z")
    axes[0].set_ylabel(r"interface $J_x$")
    axes[0].set_title("Recovered interface current")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=8)
    axes[1].plot(fluxes["z"], fluxes["left_jx"] - fluxes["right_jx"])
    axes[1].axhline(0.0, color="0.35", lw=0.8)
    axes[1].set_xlabel("z")
    axes[1].set_ylabel("left - right")
    axes[1].set_title("Interface-current jump")
    axes[1].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_fe_two_region_flux_mms_interface.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, np.ndarray]]]) -> None:
    finest = cases[-1]
    lines = [
        "# Q1 FE two-region interface-flux MMS",
        "",
        "This standalone Python reference aligns a conductivity jump with Q1",
        "element faces, solves the insulating weak form, and recovers normal",
        "current on both sides of the material interface. It does not run MUG.",
        "",
        "Potential and interface-current recovery are at roundoff because the",
        "piecewise-linear exact potential is in the conforming Q1 space and the",
        "interface is aligned to element faces.",
        "",
        "| n | max potential error | max interface target error | max interface jump |",
        "|---:|---:|---:|---:|",
    ]
    for case in cases:
        lines.append(
            f"| {int(case['x'].size)} | `{float(case['max_error']):.3e}` | "
            f"`{float(case['max_interface_target_error']):.3e}` | "
            f"`{float(case['max_interface_jump']):.3e}` |"
        )
    lines.extend(
        [
            "",
            "Finest-grid interface current:",
            f"- max interface target-current error = {float(finest['max_interface_target_error']):.3e}.",
            f"- max left/right interface-current jump = {float(finest['max_interface_jump']):.3e}.",
            "",
            "Red-team interpretation:",
            "- Interface is vertical, conforming, and exactly aligned to element faces.",
            "- Conductivity is scalar and piecewise constant.",
            "- This checks interface flux recovery in a standalone Python Q1 prototype,",
            "  not in MUG and not in an axisymmetric metric.",
        ]
    )
    (outdir / "inductionless_fe_two_region_flux_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sizes = [int(item) for item in args.sizes.split(",") if item.strip()]
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cases = [solve_case(args, size) for size in sizes]
    write_metrics(outdir, cases)
    plot_solution(outdir, cases[-1])
    plot_interface(outdir, cases[-1])
    write_summary(outdir, cases)
    print(f"Wrote Q1 FE two-region flux MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
