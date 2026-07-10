#!/usr/bin/env python3
"""Q1 finite-element inductionless potential MMS.

This standalone reference assembles the weak form

    int sigma grad(w).grad(phi) = int sigma grad(w).(u x B)

with an all-insulating natural boundary condition and a nodal mean-value gauge.
It is a sign/gauge prototype for a future electric-potential mode, not a MUG
implementation.
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
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "python"))

try:
    from OpenFUSIONToolkit.LM_MHD import (  # noqa: E402
        charge_conservation_residual_2d,
        reconstruct_inductionless_current_2d,
        wall_normal_current_extrema,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    charge_conservation_residual_2d = inductionless.charge_conservation_residual_2d
    reconstruct_inductionless_current_2d = inductionless.reconstruct_inductionless_current_2d
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--sizes", default="9,17,25,33,49", help="comma-separated grid sizes")
    parser.add_argument("--b-y", type=float, default=4.0, help="uniform transverse magnetic field")
    return parser.parse_args()


def conductivity(x: float | np.ndarray, z: float | np.ndarray) -> float | np.ndarray:
    return 1.0 + 0.25 * x + 0.15 * z


def phi_exact_raw(x: np.ndarray, z: np.ndarray) -> np.ndarray:
    return 0.02 * ((x - 0.5) ** 2 + (z - 0.5) ** 2)


def grad_phi_exact(x: float | np.ndarray, z: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    return 0.04 * (x - 0.5), 0.04 * (z - 0.5)


def target_current(x: float | np.ndarray, z: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    scale = 0.025
    return (
        scale * np.pi * np.sin(np.pi * x) * np.cos(np.pi * z),
        -scale * np.pi * np.cos(np.pi * x) * np.sin(np.pi * z),
    )


def motional_emf(x: float | np.ndarray, z: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    sigma = conductivity(x, z)
    grad_x, grad_z = grad_phi_exact(x, z)
    jx, jz = target_current(x, z)
    return grad_x + jx / sigma, grad_z + jz / sigma


def q1_shapes(xi: float, eta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shape = 0.25 * np.array(
        [
            (1.0 - xi) * (1.0 - eta),
            (1.0 + xi) * (1.0 - eta),
            (1.0 + xi) * (1.0 + eta),
            (1.0 - xi) * (1.0 + eta),
        ]
    )
    d_dxi = 0.25 * np.array([-(1.0 - eta), 1.0 - eta, 1.0 + eta, -(1.0 + eta)])
    d_deta = 0.25 * np.array([-(1.0 - xi), -(1.0 + xi), 1.0 + xi, 1.0 - xi])
    return shape, d_dxi, d_deta


def assemble_q1_system(x: np.ndarray, z: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
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
                    sigma = float(conductivity(xq, zq))
                    ex_motional, ez_motional = motional_emf(xq, zq)
                    d_dx = d_dxi * (2.0 / hx)
                    d_dz = d_deta * (2.0 / hz)
                    for i in range(4):
                        local_rhs[i] += sigma * (d_dx[i] * ex_motional + d_dz[i] * ez_motional) * det_j
                        for j in range(4):
                            local_k[i, j] += sigma * (d_dx[i] * d_dx[j] + d_dz[i] * d_dz[j]) * det_j
            for i, row in enumerate(local_nodes):
                rhs[row] += local_rhs[i]
                for j, col in enumerate(local_nodes):
                    rows.append(row)
                    cols.append(col)
                    data.append(local_k[i, j])
    return sp.csr_matrix((data, (rows, cols)), shape=(size, size)), rhs


def solve_with_nodal_mean_gauge(matrix: sp.csr_matrix, rhs: np.ndarray, mean_value: float = 0.0) -> np.ndarray:
    size = rhs.size
    constraint = sp.csr_matrix(np.ones((1, size), dtype=np.float64))
    augmented = sp.bmat([[matrix, constraint.T], [constraint, None]], format="csr")
    full_rhs = np.concatenate((rhs, np.array([mean_value * size])))
    solution = spla.spsolve(augmented, full_rhs)
    return solution[:size]


def solve_case(n: int, b_y: float) -> dict[str, np.ndarray | float | dict[str, float]]:
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.0, 1.0, n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    matrix, rhs = assemble_q1_system(x, z)
    phi = solve_with_nodal_mean_gauge(matrix, rhs).reshape(n, n)
    phi_exact = phi_exact_raw(x_grid, z_grid)
    phi_exact -= np.mean(phi_exact)
    potential_error = phi - phi_exact
    sigma_grid = conductivity(x_grid, z_grid)
    emf_x, emf_z = motional_emf(x_grid, z_grid)
    velocity = np.zeros((n, n, 3), dtype=np.float64)
    velocity[..., 0] = emf_z / b_y
    velocity[..., 2] = -emf_x / b_y
    magnetic_field = np.zeros_like(velocity)
    magnetic_field[..., 1] = b_y
    current = reconstruct_inductionless_current_2d(phi, velocity, magnetic_field, sigma_grid, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)
    boundary_flux = boundary_face_current(phi, x, z)
    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "conductivity": sigma_grid,
        "phi": phi,
        "phi_exact": phi_exact,
        "potential_error": potential_error,
        "current": current,
        "residual": residual,
        "wall_current": wall_current,
        "h": float(1.0 / (n - 1)),
        "max_error": float(np.max(np.abs(potential_error))),
        "rms_error": float(np.sqrt(np.mean(potential_error**2))),
        "max_divJ": float(np.max(np.abs(residual[2:-2, 2:-2]))),
        "max_wall_current": float(max(wall_current.values())),
        "max_fe_boundary_current": float(np.max(np.abs(boundary_flux))),
        "rms_fe_boundary_current": float(np.sqrt(np.mean(boundary_flux**2))),
    }


def fe_gradient_at(
    phi: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    ez: int,
    ex: int,
    xi: float,
    eta: float,
) -> tuple[float, float]:
    _, d_dxi, d_deta = q1_shapes(xi, eta)
    hx = x[ex + 1] - x[ex]
    hz = z[ez + 1] - z[ez]
    local_phi = np.array([phi[ez, ex], phi[ez, ex + 1], phi[ez + 1, ex + 1], phi[ez + 1, ex]])
    d_dx = d_dxi * (2.0 / hx)
    d_dz = d_deta * (2.0 / hz)
    return float(np.dot(local_phi, d_dx)), float(np.dot(local_phi, d_dz))


def boundary_face_current(phi: np.ndarray, x: np.ndarray, z: np.ndarray) -> np.ndarray:
    gauss = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))
    nx = x.size
    nz = z.size
    currents: list[float] = []
    for ez in range(nz - 1):
        z0 = z[ez]
        z1 = z[ez + 1]
        for eta in gauss:
            zq = 0.5 * (z0 + z1) + 0.5 * (z1 - z0) * eta
            grad_x, _ = fe_gradient_at(phi, x, z, ez, 0, -1.0, eta)
            emf_x, _ = motional_emf(x[0], zq)
            currents.append(float(conductivity(x[0], zq) * (grad_x - emf_x)))
            grad_x, _ = fe_gradient_at(phi, x, z, ez, nx - 2, 1.0, eta)
            emf_x, _ = motional_emf(x[-1], zq)
            currents.append(float(conductivity(x[-1], zq) * (-grad_x + emf_x)))
    for ex in range(nx - 1):
        x0 = x[ex]
        x1 = x[ex + 1]
        for xi in gauss:
            xq = 0.5 * (x0 + x1) + 0.5 * (x1 - x0) * xi
            _, grad_z = fe_gradient_at(phi, x, z, 0, ex, xi, -1.0)
            _, emf_z = motional_emf(xq, z[0])
            currents.append(float(conductivity(xq, z[0]) * (grad_z - emf_z)))
            _, grad_z = fe_gradient_at(phi, x, z, nz - 2, ex, xi, 1.0)
            _, emf_z = motional_emf(xq, z[-1])
            currents.append(float(conductivity(xq, z[-1]) * (-grad_z + emf_z)))
    return np.asarray(currents, dtype=np.float64)


def observed_orders(errors: list[float], hs: list[float]) -> list[float | None]:
    orders: list[float | None] = [None]
    for i in range(1, len(errors)):
        orders.append(np.log(errors[i - 1] / errors[i]) / np.log(hs[i - 1] / hs[i]))
    return orders


def write_metrics(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, float]]]) -> None:
    errors = [float(case["max_error"]) for case in cases]
    hs = [float(case["h"]) for case in cases]
    orders = observed_orders(errors, hs)
    with (outdir / "inductionless_fe_potential_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "n",
                "h",
                "max_error",
                "rms_error",
                "observed_order",
                "max_divJ",
                "max_wall_current",
                "max_fe_boundary_current",
                "rms_fe_boundary_current",
            ]
        )
        for case, order in zip(cases, orders):
            writer.writerow(
                [
                    int(case["x"].size),
                    f"{float(case['h']):.12e}",
                    f"{float(case['max_error']):.12e}",
                    f"{float(case['rms_error']):.12e}",
                    "" if order is None else f"{order:.6f}",
                    f"{float(case['max_divJ']):.12e}",
                    f"{float(case['max_wall_current']):.12e}",
                    f"{float(case['max_fe_boundary_current']):.12e}",
                    f"{float(case['rms_fe_boundary_current']):.12e}",
                ]
            )


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
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
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Q1 weak-form inductionless potential MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_fe_potential_mms_solution.png", dpi=180)
    plt.close(fig)


def plot_convergence(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, float]]]) -> None:
    hs = np.array([float(case["h"]) for case in cases])
    errors = np.array([float(case["max_error"]) for case in cases])
    boundary = np.array([float(case["max_fe_boundary_current"]) for case in cases])
    reference_phi = errors[-1] * (hs / hs[-1]) ** 2
    reference_flux = boundary[-1] * (hs / hs[-1])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    axes[0].loglog(hs, errors, "o-", label="Q1 FE max error")
    axes[0].loglog(hs, reference_phi, "--", label=r"$O(h^2)$ reference")
    axes[0].invert_xaxis()
    axes[0].set_xlabel("h")
    axes[0].set_ylabel(r"max $|\phi-\phi_{exact}|$")
    axes[0].set_title("Potential convergence")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=8)
    axes[1].loglog(hs, boundary, "o-", label=r"max boundary $|J\cdot n|$")
    axes[1].loglog(hs, reference_flux, "--", label=r"$O(h)$ reference")
    axes[1].invert_xaxis()
    axes[1].set_xlabel("h")
    axes[1].set_ylabel(r"boundary $|J\cdot n|$")
    axes[1].set_title("Recovered boundary flux")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_fe_potential_mms_convergence.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, float]]]) -> None:
    errors = [float(case["max_error"]) for case in cases]
    hs = [float(case["h"]) for case in cases]
    orders = observed_orders(errors, hs)
    lines = [
        "# Q1 finite-element inductionless potential MMS",
        "",
        "This standalone Python reference assembles the weak form",
        "`int sigma grad(w).grad(phi) = int sigma grad(w).(u x B)` with",
        "all-insulating natural current closure and a nodal mean-value gauge.",
        "It does not run MUG.",
        "",
        "Convergence:",
        "",
        "| n | max error | observed order |",
        "|---:|---:|---:|",
    ]
    for case, order in zip(cases, orders):
        order_text = "" if order is None else f"`{order:.3f}`"
        lines.append(f"| {int(case['x'].size)} | `{float(case['max_error']):.3e}` | {order_text} |")
    finest = cases[-1]
    lines.extend(
        [
            "",
            "Finest-grid diagnostic:",
            f"- max interior |div J| from nodal postprocessing = {float(finest['max_divJ']):.3e}.",
            f"- max wall-normal current from nodal postprocessing = {float(finest['max_wall_current']):.3e}.",
            f"- max boundary-face |J.n| from FE flux recovery = {float(finest['max_fe_boundary_current']):.3e}.",
            "- The recovered boundary-face current decreases approximately first order",
            "  over this sweep, as expected for simple Q1 element-gradient recovery.",
            "",
            "Red-team interpretation:",
            "- This checks weak-form sign, gauge, and convergence on a structured Q1 mesh.",
            "- The boundary-face flux is recovered from FE element gradients at Gauss",
            "  points, but no interface flux recovery is implemented yet.",
            "- The nodal current and wall-current diagnostics are still postprocessing",
            "  checks rather than the finite-element residual operator.",
            "- No Fortran field structure, conducting-wall Robin term, tensor conductivity,",
            "  or axisymmetric R metric is implemented here.",
        ]
    )
    (outdir / "inductionless_fe_potential_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    sizes = [int(item) for item in args.sizes.split(",") if item.strip()]
    if any(size < 3 for size in sizes):
        raise ValueError("--sizes must contain grid sizes >= 3")
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cases = [solve_case(size, args.b_y) for size in sizes]
    write_metrics(outdir, cases)
    plot_solution(outdir, cases[-1])
    plot_convergence(outdir, cases)
    write_summary(outdir, cases)
    print(f"Wrote Q1 FE inductionless potential MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
