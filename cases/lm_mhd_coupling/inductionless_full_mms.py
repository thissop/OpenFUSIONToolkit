#!/usr/bin/env python3
"""Full constant-sigma inductionless potential MMS diagnostic.

This connects the Python reference pieces for a future inductionless LM-MHD
mode: manufactured velocity, motional EMF, Neumann wall data, scalar potential
solve, reconstructed current, charge conservation, Joule heating, and Lorentz
force. It is not a MUG solve.
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
        insulating_wall_normal_gradient,
        joule_heating_density,
        lorentz_force_density,
        motional_electric_field,
        potential_source_from_motional_emf,
        reconstruct_inductionless_current_2d,
        solve_potential_neumann_2d,
        structured_gradient_2d,
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
    insulating_wall_normal_gradient = inductionless.insulating_wall_normal_gradient
    joule_heating_density = inductionless.joule_heating_density
    lorentz_force_density = inductionless.lorentz_force_density
    motional_electric_field = inductionless.motional_electric_field
    potential_source_from_motional_emf = inductionless.potential_source_from_motional_emf
    reconstruct_inductionless_current_2d = inductionless.reconstruct_inductionless_current_2d
    solve_potential_neumann_2d = inductionless.solve_potential_neumann_2d
    structured_gradient_2d = inductionless.structured_gradient_2d
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--n", type=int, default=81, help="grid points in each direction")
    parser.add_argument("--sigma", type=float, default=0.8e6, help="constant conductivity [S/m]")
    parser.add_argument("--b-y", type=float, default=5.0, help="uniform transverse field [T]")
    parser.add_argument("--j-scale", type=float, default=1.0e5, help="target current scale [A/m^2]")
    return parser.parse_args()


def build_case(args: argparse.Namespace) -> dict[str, np.ndarray | float | dict[str, float]]:
    x = np.linspace(0.0, 1.0, args.n)
    z = np.linspace(0.0, 1.0, args.n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")

    phi_exact = 0.02 * ((x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2)
    phi_exact -= np.mean(phi_exact)
    streamfunction = (args.j_scale / np.pi) * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)

    emf = structured_gradient_2d(phi_exact, x, z) + target_current / args.sigma
    magnetic_field = np.array([0.0, args.b_y, 0.0])
    velocity = np.zeros_like(target_current)
    velocity[..., 0] = emf[..., 2] / args.b_y
    velocity[..., 2] = -emf[..., 0] / args.b_y

    source = potential_source_from_motional_emf(motional_electric_field(velocity, magnetic_field), x, z)
    wall_gradient = insulating_wall_normal_gradient(emf, x, z)
    phi = solve_potential_neumann_2d(source, x, z, mean_value=0.0, normal_gradient=wall_gradient)
    current = reconstruct_inductionless_current_2d(phi, velocity, magnetic_field, args.sigma, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)
    force = lorentz_force_density(current, magnetic_field)
    joule = joule_heating_density(current, args.sigma)

    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "phi_exact": phi_exact,
        "phi": phi,
        "potential_error": phi - phi_exact,
        "target_current": target_current,
        "current": current,
        "current_error": current - target_current,
        "velocity": velocity,
        "emf": emf,
        "source": source,
        "residual": residual,
        "wall_current": wall_current,
        "force": force,
        "joule": joule,
    }


def metric_values(args: argparse.Namespace, case: dict[str, np.ndarray | float | dict[str, float]]) -> dict[str, float]:
    potential_error = case["potential_error"]
    current_error = case["current_error"]
    current = case["current"]
    velocity = case["velocity"]
    residual = case["residual"]
    wall_current = case["wall_current"]
    force = case["force"]
    joule = case["joule"]
    assert isinstance(wall_current, dict)
    return {
        "n": float(args.n),
        "sigma_S_m": args.sigma,
        "B_y_T": args.b_y,
        "max_abs_phi_error": float(np.max(np.abs(potential_error))),
        "rms_phi_error": float(np.sqrt(np.mean(potential_error**2))),
        "max_abs_current_error_A_m2": float(np.max(np.abs(current_error))),
        "max_abs_divJ_A_m3": float(np.max(np.abs(residual[2:-2, 2:-2]))),
        "rms_divJ_A_m3": float(np.sqrt(np.mean(residual[2:-2, 2:-2] ** 2))),
        "max_wall_normal_current_A_m2": float(max(wall_current.values())),
        "max_current_A_m2": float(np.max(np.linalg.norm(current, axis=-1))),
        "max_velocity_m_s": float(np.max(np.linalg.norm(velocity, axis=-1))),
        "max_lorentz_force_N_m3": float(np.max(np.linalg.norm(force, axis=-1))),
        "mean_joule_heating_W_m3": float(np.mean(joule)),
    }


def write_metrics(outdir: Path, metrics: dict[str, float]) -> None:
    with (outdir / "inductionless_full_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, f"{value:.12e}"])


def plot_potential(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharex=True, sharey=True)
    panels = [
        (case["phi_exact"], r"exact $\phi$", "viridis"),
        (case["phi"], r"solved $\phi$", "viridis"),
        (case["potential_error"], r"$\phi$ error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Full inductionless potential MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_full_mms_potential.png", dpi=180)
    plt.close(fig)


def plot_current(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    current = case["current"]
    target_current = case["target_current"]
    residual = case["residual"]
    joule = case["joule"]
    current_mag = np.linalg.norm(current, axis=-1)
    target_mag = np.linalg.norm(target_current, axis=-1)
    current_error_mag = np.linalg.norm(case["current_error"], axis=-1)

    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.2), sharex=True, sharey=True)
    panels = [
        (target_mag, r"target $|J|$ [A/m$^2$]", "viridis"),
        (current_mag, r"reconstructed $|J|$ [A/m$^2$]", "viridis"),
        (current_error_mag, r"$|J-J_{target}|$ [A/m$^2$]", "magma"),
        (residual, r"$\nabla \cdot J$ [A/m$^3$]", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes.flat, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Reconstructed inductionless current diagnostics")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_full_mms_current.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.3, 4.2))
    mid_z = current.shape[0] // 2
    ax.plot(case["x"], current[mid_z, :, 0], label=r"$J_x$")
    ax.plot(case["x"], current[mid_z, :, 2], label=r"$J_z$")
    ax.plot(case["x"], joule[mid_z, :], "--", label=r"$J^2/\sigma$")
    ax.set_xlabel("x at z=0.5")
    ax.set_title("Midplane current and Joule heating")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_full_mms_midplane.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, metrics: dict[str, float]) -> None:
    lines = [
        "# Full inductionless potential MMS",
        "",
        "This is a constant-conductivity Python reference diagnostic for the future",
        "electric-potential LM-MHD mode. It does not run MUG.",
        "",
        "The manufactured construction sets `u x B = grad(phi) + J_target/sigma`,",
        "where `J_target` is streamfunction-generated, divergence-free, and has zero",
        "normal current at the walls. The scalar solve uses",
        "`laplacian(phi) = div(u x B)` with insulating wall data",
        "`dphi/dn = (u x B).n`.",
        "",
        "Headline metrics:",
        f"- max |phi error| = {metrics['max_abs_phi_error']:.3e}.",
        f"- max |J - J_target| = {metrics['max_abs_current_error_A_m2']:.3e} A/m^2.",
        f"- max interior |div J| = {metrics['max_abs_divJ_A_m3']:.3e} A/m^3.",
        f"- max wall-normal current = {metrics['max_wall_normal_current_A_m2']:.3e} A/m^2.",
        "",
        "Generated files:",
        "- `inductionless_full_mms_metrics.csv`",
        "- `inductionless_full_mms_potential.png`",
        "- `inductionless_full_mms_current.png`",
        "- `inductionless_full_mms_midplane.png`",
        "",
        "Red-team interpretation:",
        "- Constant sigma only; no variable-coefficient or material-interface solve.",
        "- Uniform-grid finite differences only; no finite-element weak form.",
        "- No conducting-wall Robin closure.",
        "- This is a reference target for future solver work, not a solver validation.",
    ]
    (outdir / "inductionless_full_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    case = build_case(args)
    metrics = metric_values(args, case)
    write_metrics(outdir, metrics)
    plot_potential(outdir, case)
    plot_current(outdir, case)
    write_summary(outdir, metrics)
    print(f"Wrote full inductionless MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
