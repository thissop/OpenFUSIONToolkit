#!/usr/bin/env python3
"""Axisymmetric loop-voltage-driven inductionless feedback scan.

This diagnostic is the loop-voltage counterpart to
`axisymmetric_inductionless_toroidal_feedback.py`: it derives
E_phi(R) = V_loop/(2*pi*R), reconstructs inductionless current, writes a
CouplingSnapshot, and estimates toroidal-current magnetic feedback. It is not a
MUG or TokaMaker solve.
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
        CouplingSnapshot,
        axisymmetric_area_weights,
        axisymmetric_boundary_current_integrals,
        axisymmetric_cell_volumes,
        axisymmetric_charge_conservation_residual,
        axisymmetric_current_from_flux_function,
        axisymmetric_electric_field_from_potential,
        axisymmetric_integrated_power,
        axisymmetric_integrated_torque,
        axisymmetric_net_boundary_current,
        axisymmetric_wall_normal_current_extrema,
        electric_power_density,
        feedback_metrics_from_snapshot,
        inductionless_power_balance_residual,
        joule_heating_density,
        lorentz_force_density,
        mechanical_power_density,
        reconstruct_axisymmetric_inductionless_current,
        read_coupling_snapshot,
        toroidal_electric_field_from_loop_voltage,
        write_coupling_snapshot,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    coupling_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    coupling_spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(coupling_spec)
    sys.modules[coupling_spec.name] = coupling
    coupling_spec.loader.exec_module(coupling)

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    inductionless_spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(inductionless_spec)
    sys.modules[inductionless_spec.name] = inductionless
    inductionless_spec.loader.exec_module(inductionless)

    CouplingSnapshot = coupling.CouplingSnapshot
    axisymmetric_area_weights = coupling.axisymmetric_area_weights
    feedback_metrics_from_snapshot = coupling.feedback_metrics_from_snapshot
    read_coupling_snapshot = coupling.read_coupling_snapshot
    write_coupling_snapshot = coupling.write_coupling_snapshot
    axisymmetric_boundary_current_integrals = inductionless.axisymmetric_boundary_current_integrals
    axisymmetric_cell_volumes = inductionless.axisymmetric_cell_volumes
    axisymmetric_charge_conservation_residual = inductionless.axisymmetric_charge_conservation_residual
    axisymmetric_current_from_flux_function = inductionless.axisymmetric_current_from_flux_function
    axisymmetric_electric_field_from_potential = inductionless.axisymmetric_electric_field_from_potential
    axisymmetric_integrated_power = inductionless.axisymmetric_integrated_power
    axisymmetric_integrated_torque = inductionless.axisymmetric_integrated_torque
    axisymmetric_net_boundary_current = inductionless.axisymmetric_net_boundary_current
    axisymmetric_wall_normal_current_extrema = inductionless.axisymmetric_wall_normal_current_extrema
    electric_power_density = inductionless.electric_power_density
    joule_heating_density = inductionless.joule_heating_density
    inductionless_power_balance_residual = inductionless.inductionless_power_balance_residual
    lorentz_force_density = inductionless.lorentz_force_density
    mechanical_power_density = inductionless.mechanical_power_density
    reconstruct_axisymmetric_inductionless_current = inductionless.reconstruct_axisymmetric_inductionless_current
    toroidal_electric_field_from_loop_voltage = inductionless.toroidal_electric_field_from_loop_voltage


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--nr", type=int, default=81, help="radial grid points")
    parser.add_argument("--nz", type=int, default=65, help="vertical grid points")
    parser.add_argument("--sigma0", type=float, default=8.0e5, help="representative electrical conductivity [S/m]")
    parser.add_argument("--b-phi", type=float, default=5.0, help="dominant toroidal magnetic field [T]")
    parser.add_argument("--reference-field", type=float, default=5.0, help="feedback normalization field [T]")
    parser.add_argument("--selected-loop-voltage", type=float, default=0.2, help="selected map case [V]")
    parser.add_argument(
        "--loop-voltage-values",
        default="0,0.02,0.05,0.1,0.2,0.5",
        help="comma-separated loop-voltage scan values [V]",
    )
    return parser.parse_args()


def grid(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r = np.linspace(1.0, 2.0, args.nr)
    z = np.linspace(-0.5, 0.5, args.nz)
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    return r, z, r_grid, z_grid


def conductivity(args: argparse.Namespace, r_grid: np.ndarray, z_grid: np.ndarray) -> np.ndarray:
    return args.sigma0 * (1.0 + 0.08 * (r_grid - 1.5) + 0.04 * z_grid)


def probe_points() -> np.ndarray:
    return np.array(
        [[0.85, 0.0], [1.10, 0.35], [1.45, 0.0], [1.90, -0.35], [2.15, 0.0]],
        dtype=np.float64,
    )


def build_case(args: argparse.Namespace, loop_voltage: float) -> dict[str, np.ndarray | float]:
    r, z, r_grid, z_grid = grid(args)
    sigma = conductivity(args, r_grid, z_grid)
    phi = 0.02 * ((r_grid - 1.5) ** 2 + z_grid**2)
    phi -= np.mean(phi)
    flux = 2.0e4 * (r_grid - r[0]) ** 2 * (r[-1] - r_grid) ** 2 * (z_grid - z[0]) ** 2 * (z[-1] - z_grid) ** 2
    target_poloidal_current = axisymmetric_current_from_flux_function(flux, r, z)

    grad_r = np.gradient(phi, r, axis=1, edge_order=2)
    grad_z = np.gradient(phi, z, axis=0, edge_order=2)
    emf_r = grad_r + target_poloidal_current[..., 0] / sigma
    emf_z = grad_z + target_poloidal_current[..., 2] / sigma

    magnetic_field = np.zeros(r_grid.shape + (3,), dtype=np.float64)
    magnetic_field[..., 0] = 0.08 * np.sin(np.pi * (z_grid - z[0]) / (z[-1] - z[0]))
    magnetic_field[..., 1] = args.b_phi
    magnetic_field[..., 2] = 0.18 + 0.03 * (r_grid - 1.5)

    velocity = np.zeros_like(magnetic_field)
    velocity[..., 0] = emf_z / args.b_phi
    velocity[..., 2] = -emf_r / args.b_phi
    e_phi = toroidal_electric_field_from_loop_voltage(loop_voltage, r_grid)
    electric_field = axisymmetric_electric_field_from_potential(
        phi,
        r,
        z,
        toroidal_electric_field=e_phi,
    )
    current = reconstruct_axisymmetric_inductionless_current(
        phi,
        velocity,
        magnetic_field,
        sigma,
        r,
        z,
        toroidal_electric_field=e_phi,
    )
    residual = axisymmetric_charge_conservation_residual(current, r, z)
    wall_current = axisymmetric_wall_normal_current_extrema(current)
    boundary_current = axisymmetric_boundary_current_integrals(current, r, z)
    net_boundary_current = axisymmetric_net_boundary_current(current, r, z)
    force = lorentz_force_density(current, magnetic_field)
    heating = joule_heating_density(current, sigma)
    mechanical_power = mechanical_power_density(force, velocity)
    electric_power = electric_power_density(current, electric_field)
    power_balance_residual = inductionless_power_balance_residual(current, velocity, magnetic_field, sigma, electric_field)
    integrated_joule_power = axisymmetric_integrated_power(heating, r, z)
    integrated_mechanical_power = axisymmetric_integrated_power(mechanical_power, r, z)
    integrated_electric_power = axisymmetric_integrated_power(electric_power, r, z)
    integrated_power_balance_residual = axisymmetric_integrated_power(power_balance_residual, r, z)
    net_axisymmetric_torque = axisymmetric_integrated_torque(force, r, z)
    cell_volume = axisymmetric_cell_volumes(r, z)

    snapshot = CouplingSnapshot(
        points_rz=np.column_stack((r_grid.ravel(), z_grid.ravel())),
        current_density=current.reshape(-1, 3),
        time=0.0,
        cell_volume=cell_volume.ravel(),
        metadata={
            "source": "axisymmetric_loop_voltage_feedback.py",
            "loop_voltage_V": float(loop_voltage),
            "interpretation": "manufactured loop-voltage inductionless current diagnostic",
        },
    )
    return {
        "r": r,
        "z": z,
        "r_grid": r_grid,
        "z_grid": z_grid,
        "conductivity": sigma,
        "e_phi": e_phi,
        "potential": phi,
        "magnetic_field": magnetic_field,
        "velocity": velocity,
        "electric_field": electric_field,
        "current": current,
        "residual": residual,
        "wall_current": wall_current,
        "boundary_current": boundary_current,
        "net_boundary_current": net_boundary_current,
        "force": force,
        "joule_heating": heating,
        "mechanical_power": mechanical_power,
        "electric_power": electric_power,
        "power_balance_residual": power_balance_residual,
        "integrated_joule_power": integrated_joule_power,
        "integrated_mechanical_power": integrated_mechanical_power,
        "integrated_electric_power": integrated_electric_power,
        "integrated_power_balance_residual": integrated_power_balance_residual,
        "net_axisymmetric_torque": net_axisymmetric_torque,
        "cell_volume": cell_volume,
        "snapshot": snapshot,
        "loop_voltage": float(loop_voltage),
    }


def metric_row(args: argparse.Namespace, loop_voltage: float) -> dict[str, float | str]:
    case = build_case(args, loop_voltage)
    snapshot = case["snapshot"]
    area_weights = axisymmetric_area_weights(snapshot)
    metrics = feedback_metrics_from_snapshot(
        snapshot,
        probe_points(),
        reference_field=args.reference_field,
        area_weights=area_weights,
        regularization=0.03,
    )
    residual = case["residual"]
    wall_current = case["wall_current"]
    current = case["current"]
    force = case["force"]
    heating = case["joule_heating"]
    return {
        "loop_voltage": float(loop_voltage),
        "e_phi_min": float(np.min(case["e_phi"])),
        "e_phi_max": float(np.max(case["e_phi"])),
        "total_toroidal_current": metrics.total_toroidal_current,
        "absolute_toroidal_current": metrics.absolute_toroidal_current,
        "max_delta_b": metrics.max_delta_b,
        "max_feedback": metrics.max_feedback,
        "rms_feedback": metrics.rms_feedback,
        "classification": metrics.classification,
        "max_abs_j_phi": float(np.max(np.abs(current[..., 1]))),
        "max_abs_force": float(np.max(np.linalg.norm(force, axis=-1))),
        "max_joule_heating": float(np.max(heating)),
        "integrated_joule_power": float(case["integrated_joule_power"]),
        "integrated_mechanical_power": float(case["integrated_mechanical_power"]),
        "integrated_electric_power": float(case["integrated_electric_power"]),
        "integrated_power_balance_residual": float(case["integrated_power_balance_residual"]),
        "max_abs_power_balance_residual": float(np.max(np.abs(case["power_balance_residual"]))),
        "net_axisymmetric_torque": float(case["net_axisymmetric_torque"]),
        "max_abs_axisym_divJ": float(np.max(np.abs(residual[2:-2, 2:-2]))),
        "max_wall_normal_current": float(max(wall_current.values())),
        "max_abs_boundary_current": float(max(abs(value) for value in case["boundary_current"].values())),
        "net_boundary_current": float(case["net_boundary_current"]),
    }


def write_scan(outdir: Path, rows: list[dict[str, float | str]]) -> None:
    with (outdir / "axisymmetric_loop_voltage_feedback_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def save_snapshot(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    path = outdir / "axisymmetric_loop_voltage_feedback_snapshot.h5"
    write_coupling_snapshot(str(path), case["snapshot"])
    roundtrip = read_coupling_snapshot(str(path))
    if roundtrip.current_density.shape != case["snapshot"].current_density.shape:
        raise RuntimeError("CouplingSnapshot HDF5 roundtrip shape mismatch")


def plot_maps(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    r_grid = case["r_grid"]
    z_grid = case["z_grid"]
    current = case["current"]
    panels = [
        (case["e_phi"], r"$E_\phi(R)$ [V/m]", "viridis"),
        (current[..., 1], r"$J_\phi$ [A/m$^2$]", "coolwarm"),
        (np.linalg.norm(case["force"], axis=-1), r"$|J\times B|$ [N/m$^3$]", "magma"),
        (case["joule_heating"], r"$|J|^2/\sigma$ [W/m$^3$]", "magma"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(r_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle(f"Loop-voltage inductionless feedback, V_loop={float(case['loop_voltage']):.3g} V")
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_loop_voltage_feedback_maps.png", dpi=180)
    plt.close(fig)


def plot_scan(outdir: Path, rows: list[dict[str, float | str]]) -> None:
    voltage = np.array([float(row["loop_voltage"]) for row in rows])
    feedback = np.array([float(row["max_feedback"]) for row in rows])
    current = np.array([float(row["total_toroidal_current"]) for row in rows])
    force = np.array([float(row["max_abs_force"]) for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.1))
    axes[0].plot(voltage, current, "o-")
    axes[0].axhline(0.0, color="0.35", lw=0.8)
    axes[0].set_xlabel(r"$V_{loop}$ [V]")
    axes[0].set_ylabel(r"$\int J_\phi dA$ [A]")
    axes[0].set_title("Integrated toroidal current")
    axes[0].grid(True, alpha=0.25)
    axes[1].semilogy(voltage, feedback, "o-")
    axes[1].axhline(1.0e-3, color="0.35", ls="--", label="two-way threshold")
    axes[1].axhline(1.0e-2, color="0.20", ls=":", label="strong threshold")
    axes[1].set_xlabel(r"$V_{loop}$ [V]")
    axes[1].set_ylabel(r"max $|\delta B|/|B_{ref}|$")
    axes[1].set_title("Feedback triage")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=8)
    axes[2].plot(voltage, force, "o-")
    axes[2].set_xlabel(r"$V_{loop}$ [V]")
    axes[2].set_ylabel(r"max $|J\times B|$ [N/m$^3$]")
    axes[2].set_title("Lorentz-force scale")
    axes[2].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_loop_voltage_feedback_scan.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, rows: list[dict[str, float | str]], selected: dict[str, np.ndarray | float]) -> None:
    selected_row = min(rows, key=lambda row: abs(float(row["loop_voltage"]) - float(selected["loop_voltage"])))
    lines = [
        "# Axisymmetric loop-voltage feedback diagnostic",
        "",
        "This manufactured diagnostic derives `E_phi(R) = V_loop/(2*pi*R)`,",
        "reconstructs inductionless current, and estimates toroidal-current",
        "magnetic feedback. It does not run MUG or TokaMaker.",
        "",
        "Selected-case metrics:",
        f"- V_loop = {float(selected_row['loop_voltage']):.3e} V.",
        f"- E_phi range = [{float(selected_row['e_phi_min']):.3e}, {float(selected_row['e_phi_max']):.3e}] V/m.",
        f"- total toroidal current = {float(selected_row['total_toroidal_current']):.3e} A.",
        f"- max |delta B|/|B_ref| = {float(selected_row['max_feedback']):.3e}.",
        f"- feedback classification = `{selected_row['classification']}`.",
        f"- max |J_phi| = {float(selected_row['max_abs_j_phi']):.3e} A/m^2.",
        f"- max |J x B| = {float(selected_row['max_abs_force']):.3e} N/m^3.",
        f"- max Joule heating = {float(selected_row['max_joule_heating']):.3e} W/m^3.",
        f"- integrated Joule power = {float(selected_row['integrated_joule_power']):.3e} W.",
        f"- integrated mechanical power f.u = {float(selected_row['integrated_mechanical_power']):.3e} W.",
        f"- integrated electric power J.E = {float(selected_row['integrated_electric_power']):.3e} W.",
        f"- integrated power-balance residual = {float(selected_row['integrated_power_balance_residual']):.3e} W.",
        f"- max pointwise power-balance residual = {float(selected_row['max_abs_power_balance_residual']):.3e} W/m^3.",
        f"- net axisymmetric torque = {float(selected_row['net_axisymmetric_torque']):.3e} N m.",
        f"- max interior axisymmetric |div J| = {float(selected_row['max_abs_axisym_divJ']):.3e}.",
        f"- max wall-normal current = {float(selected_row['max_wall_normal_current']):.3e}.",
        f"- max integrated boundary current = {float(selected_row['max_abs_boundary_current']):.3e} A.",
        f"- net integrated boundary current = {float(selected_row['net_boundary_current']):.3e} A.",
        "",
        "Generated files:",
        "- `axisymmetric_loop_voltage_feedback_scan.csv`",
        "- `axisymmetric_loop_voltage_feedback_maps.png`",
        "- `axisymmetric_loop_voltage_feedback_scan.png`",
        "- `axisymmetric_loop_voltage_feedback_snapshot.h5`",
        "",
        "Red-team interpretation:",
        "- Manufactured structured-grid reference only.",
        "- Loop voltage sets the toroidal electric field, but the potential and",
        "  velocity fields are still prescribed.",
        "- Circular-loop feedback is a triage estimate, not a Grad-Shafranov",
        "  equilibrium update.",
    ]
    (outdir / "axisymmetric_loop_voltage_feedback_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    voltages = [float(item) for item in args.loop_voltage_values.split(",") if item.strip()]
    rows = [metric_row(args, voltage) for voltage in voltages]
    selected = build_case(args, args.selected_loop_voltage)
    write_scan(outdir, rows)
    save_snapshot(outdir, selected)
    plot_maps(outdir, selected)
    plot_scan(outdir, rows)
    write_summary(outdir, rows, selected)
    print(f"Wrote axisymmetric loop-voltage feedback diagnostics to {outdir}")


if __name__ == "__main__":
    main()
