#!/usr/bin/env python3
"""Manufactured inductionless LM-MHD diagnostic.

This script builds a divergence-free manufactured current from a streamfunction,
constructs a velocity field that reproduces it through inductionless Ohm's law
with uniform B_y and zero electric potential, and writes red-team plots/metrics.
It is a Python reference artifact for sign conventions and diagnostics, not a
MUG inductionless solver.
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
        joule_heating_density,
        lorentz_force_density,
        ohms_law_current_density,
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
    divergence_free_current_from_streamfunction = inductionless.divergence_free_current_from_streamfunction
    joule_heating_density = inductionless.joule_heating_density
    lorentz_force_density = inductionless.lorentz_force_density
    ohms_law_current_density = inductionless.ohms_law_current_density
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--nx", type=int, default=97, help="number of x grid points")
    parser.add_argument("--nz", type=int, default=97, help="number of z grid points")
    parser.add_argument("--sigma", type=float, default=0.8e6, help="electrical conductivity [S/m]")
    parser.add_argument("--b-y", type=float, default=5.0, help="uniform transverse magnetic field [T]")
    parser.add_argument("--target-j", type=float, default=1.0e5, help="target peak current density scale [A/m^2]")
    return parser.parse_args()


def build_case(args: argparse.Namespace) -> dict[str, np.ndarray | float | dict[str, float]]:
    x = np.linspace(0.0, 1.0, args.nx)
    z = np.linspace(0.0, 1.0, args.nz)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    psi_scale = args.target_j / np.pi
    streamfunction = psi_scale * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)

    velocity = np.zeros_like(target_current)
    velocity[..., 0] = target_current[..., 2] / (args.sigma * args.b_y)
    velocity[..., 2] = -target_current[..., 0] / (args.sigma * args.b_y)
    magnetic_field = np.array([0.0, args.b_y, 0.0])
    reconstructed_current = ohms_law_current_density(args.sigma, velocity, magnetic_field)
    np.testing.assert_allclose(reconstructed_current, target_current, rtol=1.0e-13, atol=1.0e-8)

    residual = charge_conservation_residual_2d(reconstructed_current, x, z)
    wall_current = wall_normal_current_extrema(reconstructed_current)
    force = lorentz_force_density(reconstructed_current, magnetic_field)
    joule = joule_heating_density(reconstructed_current, args.sigma)
    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "streamfunction": streamfunction,
        "current": reconstructed_current,
        "velocity": velocity,
        "residual": residual,
        "wall_current": wall_current,
        "force": force,
        "joule": joule,
    }


def metrics(args: argparse.Namespace, case: dict[str, np.ndarray | float | dict[str, float]]) -> dict[str, float]:
    current = case["current"]
    velocity = case["velocity"]
    residual = case["residual"]
    force = case["force"]
    joule = case["joule"]
    wall_current = case["wall_current"]
    assert isinstance(wall_current, dict)
    current_mag = np.linalg.norm(current, axis=-1)
    velocity_mag = np.linalg.norm(velocity, axis=-1)
    force_mag = np.linalg.norm(force, axis=-1)
    return {
        "sigma_S_m": args.sigma,
        "B_y_T": args.b_y,
        "target_j_A_m2": args.target_j,
        "max_current_A_m2": float(np.max(current_mag)),
        "max_velocity_m_s": float(np.max(velocity_mag)),
        "max_lorentz_force_N_m3": float(np.max(force_mag)),
        "mean_joule_heating_W_m3": float(np.mean(joule)),
        "max_abs_divJ_A_m3": float(np.max(np.abs(residual[2:-2, 2:-2]))),
        "rms_divJ_A_m3": float(np.sqrt(np.mean(residual[2:-2, 2:-2] ** 2))),
        "max_wall_normal_current_A_m2": float(max(wall_current.values())),
    }


def write_metrics(outdir: Path, values: dict[str, float]) -> None:
    path = outdir / "inductionless_mms_metrics.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in values.items():
            writer.writerow([key, f"{value:.12e}"])


def plot_fields(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    current = case["current"]
    velocity = case["velocity"]
    residual = case["residual"]
    force = case["force"]
    joule = case["joule"]
    current_mag = np.linalg.norm(current, axis=-1)
    velocity_mag = np.linalg.norm(velocity, axis=-1)
    force_x = force[..., 0]

    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.2), sharex=True, sharey=True)
    panels = [
        (current_mag, r"$|J|$ [A/m$^2$]", "viridis"),
        (velocity_mag, r"$|u|$ [m/s]", "magma"),
        (residual, r"$\nabla \cdot J$ [A/m$^3$]", "coolwarm"),
        (joule, r"$J^2/\sigma$ [W/m$^3$]", "plasma"),
    ]
    for ax, (values, label, cmap) in zip(axes.flat, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(label)
        cbar = fig.colorbar(img, ax=ax)
        cbar.set_label(label)
    fig.suptitle("Manufactured inductionless Ohm's-law diagnostic")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_mms_fields.png", dpi=180)
    plt.close(fig)

    mid_z = current.shape[0] // 2
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(case["x"], current[mid_z, :, 0], label=r"$J_x$")
    ax.plot(case["x"], current[mid_z, :, 2], label=r"$J_z$")
    ax.plot(case["x"], force_x[mid_z, :], "--", label=r"$(J \times B)_x$")
    ax.set_xlabel("x at z=0.5")
    ax.set_ylabel("manufactured fields")
    ax.set_title("Midplane sign-convention slice")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_mms_midplane.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, args: argparse.Namespace, values: dict[str, float]) -> None:
    lines = [
        "# Inductionless manufactured-solution diagnostic",
        "",
        "This is a Python reference artifact for the future electric-potential LM-MHD",
        "mode. It does not run MUG and does not validate an inductionless solver.",
        "",
        "Manufactured setup:",
        "- A streamfunction creates a rectangular-domain current with zero discrete divergence.",
        "- The wall-normal current is zero on all four rectangular boundaries.",
        "- With uniform `B_y` and zero electric potential, a velocity field is chosen so",
        "  `sigma (u x B)` reconstructs the manufactured current.",
        "",
        "Headline metrics:",
        f"- sigma = {args.sigma:.3e} S/m.",
        f"- B_y = {args.b_y:.3e} T.",
        f"- max |J| = {values['max_current_A_m2']:.3e} A/m^2.",
        f"- max |u| = {values['max_velocity_m_s']:.3e} m/s.",
        f"- max interior |div J| = {values['max_abs_divJ_A_m3']:.3e} A/m^3.",
        f"- max wall-normal current = {values['max_wall_normal_current_A_m2']:.3e} A/m^2.",
        "",
        "Generated files:",
        "- `inductionless_mms_metrics.csv`",
        "- `inductionless_mms_fields.png`",
        "- `inductionless_mms_midplane.png`",
        "",
        "Red-team interpretation:",
        "- This checks signs, algebra, and diagnostics only.",
        "- It does not solve the electric-potential Poisson equation.",
        "- It does not include material jumps, conducting walls, or finite-element operators.",
        "- The next implementation step is a real Poisson reference solve with Neumann/Robin",
        "  wall-current closure and manufactured finite-element residual checks.",
    ]
    (outdir / "inductionless_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    case = build_case(args)
    values = metrics(args, case)
    write_metrics(outdir, values)
    plot_fields(outdir, case)
    write_summary(outdir, args, values)
    print(f"Wrote inductionless MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
