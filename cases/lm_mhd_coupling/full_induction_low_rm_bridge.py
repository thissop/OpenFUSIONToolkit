#!/usr/bin/env python3
"""Reduced full-induction/inductionless bridge for a 1D channel.

This diagnostic solves the steady and transient linear full-induction equation
for a prescribed parabolic channel flow in a transverse applied field and
compares the resulting Ampere current to the inductionless Ohm's-law current.
It is a reference/asymptotic bridge, not a MUG run.
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
    from OpenFUSIONToolkit.LM_MHD import MU0  # noqa: E402
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    coupling_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    coupling_spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(coupling_spec)
    sys.modules[coupling_spec.name] = coupling
    coupling_spec.loader.exec_module(coupling)
    MU0 = coupling.MU0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--half-width", type=float, default=0.05, help="channel half-width a [m]")
    parser.add_argument("--b0", type=float, default=5.0, help="transverse applied field [T]")
    parser.add_argument("--nz", type=int, default=401, help="1D grid points")
    parser.add_argument("--selected-velocity", type=float, default=0.2, help="selected profile velocity scale [m/s]")
    parser.add_argument("--selected-sigma", type=float, default=8.0e5, help="selected conductivity [S/m]")
    parser.add_argument("--velocities", default="0.02,0.05,0.1,0.2,0.5,1.0", help="velocity scan [m/s]")
    parser.add_argument("--sigmas", default="2e5,8e5,3.2e6", help="conductivity scan [S/m]")
    return parser.parse_args()


def velocity_profile(z: np.ndarray, velocity_scale: float, half_width: float) -> np.ndarray:
    return velocity_scale * (1.0 - (z / half_width) ** 2)


def bridge_solution(
    z: np.ndarray,
    velocity_scale: float,
    magnetic_field: float,
    conductivity: float,
    half_width: float,
) -> dict[str, np.ndarray | float]:
    velocity = velocity_profile(z, velocity_scale, half_width)
    mean_velocity = float(np.trapezoid(velocity, z) / (2.0 * half_width))
    electric_y = magnetic_field * mean_velocity
    inductionless_current = conductivity * (electric_y - velocity * magnetic_field)
    induced_bx = np.zeros_like(z)
    dz = np.diff(z)
    increments = 0.5 * MU0 * (inductionless_current[:-1] + inductionless_current[1:]) * dz
    induced_bx[1:] = np.cumsum(increments)
    ampere_current = np.gradient(induced_bx, z, edge_order=2) / MU0
    current_scale = max(np.finfo(float).tiny, float(np.max(np.abs(inductionless_current))))
    rm = MU0 * conductivity * velocity_scale * half_width
    induced_ratio = float(np.max(np.abs(induced_bx)) / abs(magnetic_field))
    current_error = float(np.max(np.abs(ampere_current[2:-2] - inductionless_current[2:-2])) / current_scale)
    return {
        "z": z,
        "velocity": velocity,
        "mean_velocity": mean_velocity,
        "electric_y": electric_y,
        "inductionless_current": inductionless_current,
        "induced_bx": induced_bx,
        "ampere_current": ampere_current,
        "rm": float(rm),
        "induced_ratio": induced_ratio,
        "current_error": current_error,
        "wall_bx_mismatch": float(max(abs(induced_bx[0]), abs(induced_bx[-1]))),
        "net_current": float(np.trapezoid(inductionless_current, z)),
        "magnetic_diffusion_time": float(MU0 * conductivity * half_width**2),
        "flow_time": float(half_width / velocity_scale),
    }


def transient_relaxation(
    z: np.ndarray,
    steady_bx: np.ndarray,
    velocity_scale: float,
    magnetic_field: float,
    conductivity: float,
    half_width: float,
    n_steps: int = 300,
    t_end_factor: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    eta_m = 1.0 / (MU0 * conductivity)
    tau_m = MU0 * conductivity * half_width**2
    times = np.linspace(0.0, t_end_factor * tau_m, n_steps + 1)
    dt = times[1] - times[0]
    dz = z[1] - z[0]
    interior = slice(1, -1)
    n_interior = z.size - 2
    main = np.full(n_interior, -2.0 / dz**2)
    off = np.full(n_interior - 1, 1.0 / dz**2)
    laplacian = sp.diags((off, main, off), (-1, 0, 1), format="csc")
    matrix = sp.eye(n_interior, format="csc") - dt * eta_m * laplacian
    solve = spla.factorized(matrix)
    velocity = velocity_profile(z, velocity_scale, half_width)
    forcing = magnetic_field * np.gradient(velocity, z, edge_order=2)
    bx = np.zeros(n_interior, dtype=np.float64)
    steady = steady_bx[interior]
    steady_norm = max(np.finfo(float).tiny, float(np.sqrt(np.mean(steady**2))))
    errors = [1.0]
    for _ in range(n_steps):
        bx = solve(bx + dt * forcing[interior])
        errors.append(float(np.sqrt(np.mean((bx - steady) ** 2)) / steady_norm))
    return times / tau_m, np.asarray(errors)


def scan(args: argparse.Namespace) -> tuple[np.ndarray, list[dict[str, float]]]:
    z = np.linspace(-args.half_width, args.half_width, args.nz)
    velocities = [float(item) for item in args.velocities.split(",") if item.strip()]
    sigmas = [float(item) for item in args.sigmas.split(",") if item.strip()]
    rows: list[dict[str, float]] = []
    for sigma in sigmas:
        for velocity in velocities:
            result = bridge_solution(z, velocity, args.b0, sigma, args.half_width)
            rows.append(
                {
                    "sigma": float(sigma),
                    "velocity": float(velocity),
                    "rm": float(result["rm"]),
                    "max_induced_b_over_b0": float(result["induced_ratio"]),
                    "current_linf_rel_error": float(result["current_error"]),
                    "wall_bx_mismatch": float(result["wall_bx_mismatch"]),
                    "net_current": float(result["net_current"]),
                    "magnetic_diffusion_time": float(result["magnetic_diffusion_time"]),
                    "flow_time": float(result["flow_time"]),
                    "tau_m_over_tau_u": float(result["magnetic_diffusion_time"] / result["flow_time"]),
                }
            )
    return z, rows


def write_csv(outdir: Path, rows: list[dict[str, float]]) -> None:
    with (outdir / "full_induction_low_rm_bridge_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_profiles(outdir: Path, selected: dict[str, np.ndarray | float], magnetic_field: float) -> None:
    z = selected["z"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    axes[0].plot(selected["velocity"], z)
    axes[0].axvline(selected["mean_velocity"], color="0.35", ls="--", label="mean")
    axes[0].set_xlabel("u_x [m/s]")
    axes[0].set_ylabel("z [m]")
    axes[0].set_title("Prescribed flow")
    axes[0].legend(fontsize=8)
    axes[1].plot(selected["inductionless_current"], z, label="Ohm")
    axes[1].plot(selected["ampere_current"], z, "--", label=r"$(1/\mu_0) db_x/dz$")
    axes[1].set_xlabel(r"$J_y$ [A/m$^2$]")
    axes[1].set_title("Current bridge")
    axes[1].legend(fontsize=8)
    axes[2].plot(np.asarray(selected["induced_bx"]) / magnetic_field, z)
    axes[2].set_xlabel(r"$b_x/B_0$")
    axes[2].set_title("Induced field")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "full_induction_low_rm_bridge_profiles.png", dpi=180)
    plt.close(fig)


def plot_scaling(outdir: Path, rows: list[dict[str, float]]) -> None:
    sigmas = sorted({row["sigma"] for row in rows})
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2))
    for sigma in sigmas:
        subset = sorted([row for row in rows if row["sigma"] == sigma], key=lambda row: row["rm"])
        rm = np.array([row["rm"] for row in subset])
        induced = np.array([row["max_induced_b_over_b0"] for row in subset])
        error = np.array([row["current_linf_rel_error"] for row in subset])
        axes[0].loglog(rm, induced, "o-", label=fr"$\sigma={sigma:.1e}$")
        axes[1].loglog(rm, error, "o-", label=fr"$\sigma={sigma:.1e}$")
    axes[0].loglog([1.0e-4, 1.0], [1.0e-4, 1.0], color="0.35", ls="--", label="slope 1")
    axes[0].axvline(0.1, color="0.25", ls=":", label="Rm=0.1")
    axes[0].set_xlabel(r"$R_m$")
    axes[0].set_ylabel(r"max $|b_x|/B_0$")
    axes[0].set_title("Induced-field scale")
    axes[1].axvline(0.1, color="0.25", ls=":")
    axes[1].set_xlabel(r"$R_m$")
    axes[1].set_ylabel(r"$J$ bridge rel. error")
    axes[1].set_title("Ampere vs Ohm current")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "full_induction_low_rm_bridge_scaling.png", dpi=180)
    plt.close(fig)


def plot_transient(outdir: Path, time_over_tau: np.ndarray, error: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogy(time_over_tau, error)
    ax.set_xlabel(r"$t/\tau_m$")
    ax.set_ylabel(r"$||b-b_{steady}||_2/||b_{steady}||_2$")
    ax.set_title("Magnetic diffusion relaxation")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "full_induction_low_rm_bridge_transient.png", dpi=180)
    plt.close(fig)


def write_summary(
    outdir: Path,
    args: argparse.Namespace,
    selected: dict[str, np.ndarray | float],
    rows: list[dict[str, float]],
    final_transient_error: float,
) -> None:
    row = min(
        rows,
        key=lambda item: abs(item["velocity"] - args.selected_velocity) + abs(item["sigma"] - args.selected_sigma) / args.selected_sigma,
    )
    max_rm = max(row_data["rm"] for row_data in rows)
    max_ratio = max(row_data["max_induced_b_over_b0"] for row_data in rows)
    lines = [
        "# Full-induction / inductionless low-Rm bridge",
        "",
        "This reduced 1D channel diagnostic solves the linear full-induction",
        "steady relation for induced `b_x` and compares Ampere current",
        "`(1/mu0) db_x/dz` to the inductionless Ohm's-law current with net",
        "streamwise current constrained to zero. It is not a MUG run.",
        "",
        "Selected case:",
        f"- half-width `a = {args.half_width:.3e} m`.",
        f"- B0 = {args.b0:.3e} T.",
        f"- U = {args.selected_velocity:.3e} m/s.",
        f"- sigma = {args.selected_sigma:.3e} S/m.",
        f"- Rm = {float(selected['rm']):.3e}.",
        f"- max |b_x|/B0 = {float(selected['induced_ratio']):.3e}.",
        f"- current bridge Linf relative error = {float(selected['current_error']):.3e}.",
        f"- net current integral = {float(selected['net_current']):.3e} A/m.",
        f"- wall b_x mismatch = {float(selected['wall_bx_mismatch']):.3e} T.",
        f"- final transient relative error at 3 tau_m = {final_transient_error:.3e}.",
        "",
        "Scan extrema:",
        f"- max Rm in scan = {max_rm:.3e}.",
        f"- max induced-field ratio in scan = {max_ratio:.3e}.",
        "",
        "Generated files:",
        "- `full_induction_low_rm_bridge_scan.csv`",
        "- `full_induction_low_rm_bridge_profiles.png`",
        "- `full_induction_low_rm_bridge_scaling.png`",
        "- `full_induction_low_rm_bridge_transient.png`",
        "",
        "Red-team interpretation:",
        "- Full induction is the parent model; inductionless is justified only",
        "  when induced-field ratios and Rm are small for the target regime.",
        "- This bridge uses a prescribed 1D flow and linear full-induction",
        "  equation, so it does not validate MUG or nonlinear 3D induction.",
        "- The result is a sign/asymptotic consistency check: in this controlled",
        "  low-Rm setting, Ampere current from induced `b` and Ohm's-law current",
        "  must agree.",
    ]
    (outdir / "full_induction_low_rm_bridge_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    z, rows = scan(args)
    selected = bridge_solution(z, args.selected_velocity, args.b0, args.selected_sigma, args.half_width)
    time_over_tau, transient_error = transient_relaxation(
        z,
        selected["induced_bx"],
        args.selected_velocity,
        args.b0,
        args.selected_sigma,
        args.half_width,
    )
    write_csv(outdir, rows)
    plot_profiles(outdir, selected, args.b0)
    plot_scaling(outdir, rows)
    plot_transient(outdir, time_over_tau, transient_error)
    write_summary(outdir, args, selected, rows, float(transient_error[-1]))
    print(f"Wrote full-induction low-Rm bridge diagnostics to {outdir}")


if __name__ == "__main__":
    main()
