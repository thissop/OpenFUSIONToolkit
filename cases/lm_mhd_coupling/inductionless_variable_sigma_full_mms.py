#!/usr/bin/env python3
"""Full smooth variable-conductivity inductionless MMS diagnostic.

This connects variable conductivity, motional EMF source generation, Neumann
wall flux, potential solve, current reconstruction, and charge-conservation
diagnostics. It is a Python reference artifact, not a MUG solve.
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "python"))

try:
    from OpenFUSIONToolkit.LM_MHD import (  # noqa: E402
        charge_conservation_residual_2d,
        divergence_free_current_from_streamfunction,
        insulating_wall_normal_flux,
        motional_electric_field,
        reconstruct_inductionless_current_2d,
        solve_variable_conductivity_neumann_2d,
        structured_gradient_2d,
        wall_normal_current_extrema,
        weighted_potential_source_from_motional_emf,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    charge_conservation_residual_2d = inductionless.charge_conservation_residual_2d
    divergence_free_current_from_streamfunction = inductionless.divergence_free_current_from_streamfunction
    insulating_wall_normal_flux = inductionless.insulating_wall_normal_flux
    motional_electric_field = inductionless.motional_electric_field
    reconstruct_inductionless_current_2d = inductionless.reconstruct_inductionless_current_2d
    solve_variable_conductivity_neumann_2d = inductionless.solve_variable_conductivity_neumann_2d
    structured_gradient_2d = inductionless.structured_gradient_2d
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema
    weighted_potential_source_from_motional_emf = inductionless.weighted_potential_source_from_motional_emf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--n", type=int, default=81, help="grid points in each direction")
    parser.add_argument("--b-y", type=float, default=5.0, help="uniform transverse magnetic field [T]")
    return parser.parse_args()


def build_case(args: argparse.Namespace) -> dict[str, np.ndarray | float | dict[str, float]]:
    x = np.linspace(0.0, 1.0, args.n)
    z = np.linspace(0.0, 1.0, args.n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    conductivity = 1.0 + 0.3 * x_grid + 0.2 * z_grid
    phi_exact = 0.02 * ((x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2)
    phi_exact -= np.mean(phi_exact)
    streamfunction = 0.03 * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)

    emf = structured_gradient_2d(phi_exact, x, z) + target_current / conductivity[..., None]
    magnetic_field = np.array([0.0, args.b_y, 0.0])
    velocity = np.zeros_like(target_current)
    velocity[..., 0] = emf[..., 2] / args.b_y
    velocity[..., 2] = -emf[..., 0] / args.b_y

    motional = motional_electric_field(velocity, magnetic_field)
    source = weighted_potential_source_from_motional_emf(conductivity, motional, x, z)
    normal_flux = insulating_wall_normal_flux(conductivity, motional, x, z)
    phi = solve_variable_conductivity_neumann_2d(
        conductivity,
        source,
        x,
        z,
        mean_value=0.0,
        normal_flux=normal_flux,
    )
    current = reconstruct_inductionless_current_2d(phi, velocity, magnetic_field, conductivity, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)
    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "conductivity": conductivity,
        "phi_exact": phi_exact,
        "phi": phi,
        "potential_error": phi - phi_exact,
        "target_current": target_current,
        "current": current,
        "current_error": current - target_current,
        "velocity": velocity,
        "source": source,
        "residual": residual,
        "wall_current": wall_current,
    }


def metric_values(case: dict[str, np.ndarray | float | dict[str, float]]) -> dict[str, float]:
    residual = case["residual"]
    wall_current = case["wall_current"]
    assert isinstance(wall_current, dict)
    return {
        "max_abs_phi_error": float(np.max(np.abs(case["potential_error"]))),
        "rms_phi_error": float(np.sqrt(np.mean(case["potential_error"] ** 2))),
        "max_abs_current_error": float(np.max(np.abs(case["current_error"]))),
        "rms_current_error": float(np.sqrt(np.mean(case["current_error"] ** 2))),
        "max_abs_divJ": float(np.max(np.abs(residual[2:-2, 2:-2]))),
        "rms_divJ": float(np.sqrt(np.mean(residual[2:-2, 2:-2] ** 2))),
        "max_wall_normal_current": float(max(wall_current.values())),
    }


def write_metrics(outdir: Path, metrics: dict[str, float]) -> None:
    with (outdir / "inductionless_variable_sigma_full_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, f"{value:.12e}"])


def plot_potential(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    panels = [
        (case["conductivity"], r"$\sigma(x,z)$", "viridis"),
        (case["phi_exact"], r"exact $\phi$", "magma"),
        (case["phi"], r"solved $\phi$", "magma"),
        (case["potential_error"], r"$\phi$ error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Full variable-conductivity inductionless MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_variable_sigma_full_mms_potential.png", dpi=180)
    plt.close(fig)


def plot_current(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    current = case["current"]
    target = case["target_current"]
    residual = case["residual"]
    current_error_mag = np.linalg.norm(case["current_error"], axis=-1)
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.2), sharex=True, sharey=True)
    panels = [
        (np.linalg.norm(target, axis=-1), r"target $|J|$", "viridis"),
        (np.linalg.norm(current, axis=-1), r"reconstructed $|J|$", "viridis"),
        (current_error_mag, r"$|J-J_{target}|$", "magma"),
        (residual, r"$\nabla\cdot J$", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes.flat, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Full variable-conductivity current diagnostics")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_variable_sigma_full_mms_current.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, metrics: dict[str, float]) -> None:
    lines = [
        "# Full variable-conductivity inductionless MMS",
        "",
        "This is a Python reference diagnostic connecting variable conductivity,",
        "motional EMF source generation, insulating wall flux, potential solve, and",
        "current reconstruction. It does not run MUG.",
        "",
        "Headline metrics:",
        f"- max |phi error| = {metrics['max_abs_phi_error']:.3e}.",
        f"- max |J - J_target| = {metrics['max_abs_current_error']:.3e}.",
        f"- max interior |div J| = {metrics['max_abs_divJ']:.3e}.",
        f"- max wall-normal current = {metrics['max_wall_normal_current']:.3e}.",
        "",
        "Generated files:",
        "- `inductionless_variable_sigma_full_mms_metrics.csv`",
        "- `inductionless_variable_sigma_full_mms_potential.png`",
        "- `inductionless_variable_sigma_full_mms_current.png`",
        "",
        "Red-team interpretation:",
        "- Smooth scalar conductivity only.",
        "- Uniform-grid finite differences only.",
        "- No material jump in the full inductionless reconstruction.",
        "- No conducting-wall Robin closure.",
    ]
    (outdir / "inductionless_variable_sigma_full_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    case = build_case(args)
    metrics = metric_values(case)
    write_metrics(outdir, metrics)
    plot_potential(outdir, case)
    plot_current(outdir, case)
    write_summary(outdir, metrics)
    print(f"Wrote full variable-conductivity inductionless MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
