#!/usr/bin/env python3
"""Two-region full inductionless MMS diagnostic.

This combines a face-aligned conductivity jump with motional-EMF source
generation, insulating wall flux, Neumann potential recovery, current
reconstruction, and conservative interface-current checks. It is a Python
finite-difference reference artifact, not a MUG solve.
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
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema
    weighted_potential_source_from_motional_emf = inductionless.weighted_potential_source_from_motional_emf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--nx", type=int, default=64, help="even number of x grid points")
    parser.add_argument("--nz", type=int, default=65, help="number of z grid points")
    parser.add_argument("--sigma-left", type=float, default=1.0, help="left-region conductivity")
    parser.add_argument("--sigma-right", type=float, default=4.0, help="right-region conductivity")
    parser.add_argument("--flux", type=float, default=0.7, help="constant sigma grad(phi) flux")
    parser.add_argument("--current-scale", type=float, default=0.04, help="streamfunction current amplitude")
    parser.add_argument("--b-y", type=float, default=3.0, help="uniform transverse magnetic field")
    return parser.parse_args()


def build_case(args: argparse.Namespace) -> dict[str, np.ndarray | float | dict[str, float]]:
    if args.nx % 2 != 0:
        raise ValueError("--nx must be even so the material interface lies on a grid face")
    if args.sigma_left <= 0.0 or args.sigma_right <= 0.0:
        raise ValueError("conductivities must be positive")
    if args.b_y == 0.0:
        raise ValueError("--b-y must be nonzero")

    x = np.linspace(0.0, 1.0, args.nx)
    z = np.linspace(0.0, 1.0, args.nz)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    interface = 0.5
    conductivity = np.where(x_grid < interface, args.sigma_left, args.sigma_right)
    phi_1d = np.where(
        x <= interface,
        args.flux * x / args.sigma_left,
        args.flux * interface / args.sigma_left
        + args.flux * (x - interface) / args.sigma_right,
    )
    phi_exact = np.broadcast_to(phi_1d, x_grid.shape).copy()
    phi_exact -= np.mean(phi_exact)
    streamfunction = args.current_scale * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)

    emf = np.zeros_like(target_current)
    emf[..., 0] = args.flux / conductivity + target_current[..., 0] / conductivity
    emf[..., 2] = target_current[..., 2] / conductivity
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
    left_idx = args.nx // 2 - 1
    right_idx = args.nx // 2
    mask = bulk_mask(args.nx, args.nz, left_idx, right_idx)
    interface_current, target_interface_current = interface_face_current(
        phi,
        conductivity,
        motional,
        target_current,
        x,
        left_idx,
        right_idx,
        args.sigma_left,
        args.sigma_right,
    )

    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "interface": interface,
        "left_idx": float(left_idx),
        "right_idx": float(right_idx),
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
        "bulk_mask": mask,
        "interface_current": interface_current,
        "target_interface_current": target_interface_current,
        "wall_current": wall_current,
        "sigma_left": float(args.sigma_left),
        "sigma_right": float(args.sigma_right),
        "flux": float(args.flux),
    }


def bulk_mask(nx: int, nz: int, left_idx: int, right_idx: int) -> np.ndarray:
    mask = np.ones((nz, nx), dtype=bool)
    mask[:, left_idx - 1 : right_idx + 2] = False
    mask[:2, :] = False
    mask[-2:, :] = False
    mask[:, :2] = False
    mask[:, -2:] = False
    return mask


def interface_face_current(
    phi: np.ndarray,
    conductivity: np.ndarray,
    motional: np.ndarray,
    target_current: np.ndarray,
    x: np.ndarray,
    left_idx: int,
    right_idx: int,
    sigma_left: float,
    sigma_right: float,
) -> tuple[np.ndarray, np.ndarray]:
    h = float(x[1] - x[0])
    sigma_face = 2.0 * sigma_left * sigma_right / (sigma_left + sigma_right)
    weighted_emf_x = conductivity * motional[..., 0]
    weighted_emf_face = 0.5 * (weighted_emf_x[:, left_idx] + weighted_emf_x[:, right_idx])
    face_current = -sigma_face * (phi[:, right_idx] - phi[:, left_idx]) / h + weighted_emf_face
    target_face_current = 0.5 * (target_current[:, left_idx, 0] + target_current[:, right_idx, 0])
    return face_current, target_face_current


def metric_values(case: dict[str, np.ndarray | float | dict[str, float]]) -> dict[str, float]:
    mask = case["bulk_mask"]
    wall_current = case["wall_current"]
    assert isinstance(wall_current, dict)
    current_error = case["current_error"]
    current_error_mag = np.linalg.norm(current_error, axis=-1)
    residual = case["residual"]
    interface_error = case["interface_current"] - case["target_interface_current"]
    interface_cols = ~mask
    return {
        "max_abs_phi_error": float(np.max(np.abs(case["potential_error"]))),
        "rms_phi_error": float(np.sqrt(np.mean(case["potential_error"] ** 2))),
        "max_bulk_abs_current_error": float(np.max(np.abs(current_error[mask]))),
        "rms_bulk_current_error": float(np.sqrt(np.mean(current_error[mask] ** 2))),
        "max_all_pointwise_current_error": float(np.max(current_error_mag)),
        "max_excluded_interface_band_current_error": float(np.max(current_error_mag[interface_cols])),
        "max_bulk_abs_divJ": float(np.max(np.abs(residual[mask]))),
        "max_wall_normal_current": float(max(wall_current.values())),
        "max_interface_face_current_error": float(np.max(np.abs(interface_error))),
        "rms_interface_face_current_error": float(np.sqrt(np.mean(interface_error**2))),
    }


def write_metrics(outdir: Path, metrics: dict[str, float]) -> None:
    with (outdir / "inductionless_two_region_full_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, f"{value:.12e}"])


def plot_potential(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    interface = float(case["interface"])
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    panels = [
        (case["conductivity"], r"$\sigma$", "viridis"),
        (case["phi"], r"solved $\phi$", "magma"),
        (case["potential_error"], r"$\phi$ error", "coolwarm"),
        (case["source"], r"$\nabla\cdot(\sigma u\times B)$", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.axvline(interface, color="white", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Two-region full inductionless MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_two_region_full_mms_potential.png", dpi=180)
    plt.close(fig)


def plot_current(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    interface = float(case["interface"])
    current = case["current"]
    target = case["target_current"]
    mask = case["bulk_mask"]
    residual = np.array(case["residual"], copy=True)
    current_error = np.linalg.norm(case["current_error"], axis=-1)
    bulk_error = current_error.copy()
    bulk_error[~mask] = np.nan
    residual_plot = residual.copy()
    residual_plot[~mask] = np.nan
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 7.2), sharex=True, sharey=True)
    panels = [
        (np.linalg.norm(target, axis=-1), r"target $|J|$", "viridis"),
        (np.linalg.norm(current, axis=-1), r"reconstructed $|J|$", "viridis"),
        (bulk_error, r"bulk $|J-J_{target}|$", "magma"),
        (residual_plot, r"bulk $\nabla\cdot J$", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes.flat, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.axvline(interface, color="white", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Two-region full inductionless current diagnostics")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_two_region_full_mms_current.png", dpi=180)
    plt.close(fig)


def plot_interface(outdir: Path, case: dict[str, np.ndarray | float | dict[str, float]]) -> None:
    z = case["z"]
    interface_current = case["interface_current"]
    target = case["target_interface_current"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), sharex=True)
    axes[0].plot(z, target, label="target face current")
    axes[0].plot(z, interface_current, "--", label="reconstructed face current")
    axes[0].set_xlabel("z")
    axes[0].set_ylabel(r"$J_x$ at material face")
    axes[0].set_title("Conservative interface current")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=8)

    error = interface_current - target
    axes[1].plot(z, error)
    axes[1].axhline(0.0, color="0.35", lw=0.8)
    axes[1].set_xlabel("z")
    axes[1].set_ylabel("error")
    axes[1].set_title("Face-current error")
    axes[1].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_two_region_full_mms_interface.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, metrics: dict[str, float]) -> None:
    lines = [
        "# Two-region full inductionless MMS",
        "",
        "This is a Python finite-difference reference diagnostic connecting a",
        "face-aligned conductivity jump, motional-EMF source generation,",
        "insulating wall flux, Neumann potential recovery, reconstructed current,",
        "and conservative interface-current checks. It does not run MUG.",
        "",
        "Headline metrics:",
        f"- max |phi error| = {metrics['max_abs_phi_error']:.3e}.",
        f"- max bulk |J - J_target| = {metrics['max_bulk_abs_current_error']:.3e}.",
        f"- max bulk |div J| = {metrics['max_bulk_abs_divJ']:.3e}.",
        f"- max wall-normal current = {metrics['max_wall_normal_current']:.3e}.",
        f"- max conservative interface-face current error = {metrics['max_interface_face_current_error']:.3e}.",
        f"- max all-point pointwise current error = {metrics['max_all_pointwise_current_error']:.3e}.",
        "",
        "Generated files:",
        "- `inductionless_two_region_full_mms_metrics.csv`",
        "- `inductionless_two_region_full_mms_potential.png`",
        "- `inductionless_two_region_full_mms_current.png`",
        "- `inductionless_two_region_full_mms_interface.png`",
        "",
        "Red-team interpretation:",
        "- The interface is face-aligned and exactly vertical.",
        "- Bulk pointwise current errors exclude a narrow interface band because",
        "  nodal central differences are not a conservative current operator at",
        "  discontinuous conductivity jumps.",
        "- The interface claim is therefore based on a conservative face-current",
        "  reconstruction using the same harmonic-face conductivity as the",
        "  reference potential operator.",
        "- No tensor conductivity, curved material boundary, conducting-wall",
        "  Robin closure, finite-element weak form, or MUG solve is included.",
    ]
    (outdir / "inductionless_two_region_full_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    case = build_case(args)
    metrics = metric_values(case)
    write_metrics(outdir, metrics)
    plot_potential(outdir, case)
    plot_current(outdir, case)
    plot_interface(outdir, case)
    write_summary(outdir, metrics)
    print(f"Wrote two-region full inductionless MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
