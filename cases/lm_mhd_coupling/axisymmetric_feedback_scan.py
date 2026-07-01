#!/usr/bin/env python3
"""Axisymmetric toroidal-current feedback diagnostic for LM-MHD coupling.

This script does not run MUG or TokaMaker. It creates a manufactured liquid-
metal toroidal-current patch, round-trips it through the CouplingSnapshot HDF5
schema, and estimates the poloidal-field perturbation at plasma-side probes
using circular-loop Biot-Savart geometry. The result is a triage artifact for
deciding whether two-way plasma/liquid-metal coupling deserves a real coupled
solve.
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
        axisymmetric_toroidal_current_field,
        feedback_metrics_from_snapshot,
        read_coupling_snapshot,
        write_coupling_snapshot,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    coupling_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = coupling
    spec.loader.exec_module(coupling)
    CouplingSnapshot = coupling.CouplingSnapshot
    axisymmetric_area_weights = coupling.axisymmetric_area_weights
    axisymmetric_toroidal_current_field = coupling.axisymmetric_toroidal_current_field
    feedback_metrics_from_snapshot = coupling.feedback_metrics_from_snapshot
    read_coupling_snapshot = coupling.read_coupling_snapshot
    write_coupling_snapshot = coupling.write_coupling_snapshot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--b-ref", type=float, default=5.0, help="reference field for feedback ratio [T]")
    parser.add_argument(
        "--baseline-j",
        type=float,
        default=1.0e5,
        help="baseline peak toroidal current density for map plots [A/m^2]",
    )
    return parser.parse_args()


def blanket_grid(nr: int = 62, nz: int = 74) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r_edges = np.linspace(1.45, 2.15, nr + 1)
    z_edges = np.linspace(-0.55, 0.55, nz + 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    r_grid, z_grid = np.meshgrid(r_centers, z_centers, indexing="xy")
    area = np.full(r_grid.size, (r_edges[1] - r_edges[0]) * (z_edges[1] - z_edges[0]))
    return r_grid, z_grid, area, np.column_stack((r_grid.ravel(), z_grid.ravel()))


def normalized_current_shape(r_grid: np.ndarray, z_grid: np.ndarray) -> np.ndarray:
    core = np.exp(-((r_grid - 1.75) / 0.22) ** 2 - (z_grid / 0.30) ** 2)
    upper_edge = np.exp(-((r_grid - 2.02) / 0.08) ** 2 - ((z_grid - 0.36) / 0.10) ** 2)
    lower_edge = np.exp(-((r_grid - 2.02) / 0.08) ** 2 - ((z_grid + 0.36) / 0.10) ** 2)
    shape = core + 0.35 * (upper_edge + lower_edge)
    return shape / np.max(shape)


def make_snapshot(peak_j_phi: float, r_grid: np.ndarray, z_grid: np.ndarray, area: np.ndarray) -> CouplingSnapshot:
    shape = normalized_current_shape(r_grid, z_grid)
    points = np.column_stack((r_grid.ravel(), z_grid.ravel()))
    current_density = np.zeros((points.shape[0], 3))
    current_density[:, 1] = peak_j_phi * shape.ravel()
    volume = 2.0 * np.pi * points[:, 0] * area
    return CouplingSnapshot(
        points_rz=points,
        current_density=current_density,
        time=0.0,
        cell_volume=volume,
        metadata={
            "source": "manufactured_axisymmetric_feedback_scan",
            "peak_j_phi_A_m2": float(peak_j_phi),
            "scope": "triage_not_solver_validation",
        },
    )


def probe_grid(nr: int = 70, nz: int = 80) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.linspace(0.75, 1.30, nr)
    z = np.linspace(-0.65, 0.65, nz)
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    probes = np.column_stack((r_grid.ravel(), z_grid.ravel()))
    return r_grid, z_grid, probes


def run_scan(args: argparse.Namespace, outdir: Path) -> list[dict[str, float | str]]:
    r_grid, z_grid, area, _ = blanket_grid()
    _, _, probes = probe_grid()
    amplitudes = np.array([1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6, 3.0e6])
    rows: list[dict[str, float | str]] = []
    for amplitude in amplitudes:
        snap = make_snapshot(amplitude, r_grid, z_grid, area)
        metrics = feedback_metrics_from_snapshot(snap, probes, reference_field=args.b_ref)
        rows.append(
            {
                "peak_j_phi_A_m2": float(amplitude),
                "total_toroidal_current_A": metrics.total_toroidal_current,
                "absolute_toroidal_current_A": metrics.absolute_toroidal_current,
                "centroid_R_m": metrics.current_centroid_rz[0],
                "centroid_Z_m": metrics.current_centroid_rz[1],
                "max_delta_B_T": metrics.max_delta_b,
                "rms_delta_B_T": metrics.rms_delta_b,
                "max_feedback": metrics.max_feedback,
                "rms_feedback": metrics.rms_feedback,
                "classification": metrics.classification,
            }
        )

    baseline_snapshot = make_snapshot(args.baseline_j, r_grid, z_grid, area)
    h5_path = outdir / "axisymmetric_feedback_snapshot.h5"
    write_coupling_snapshot(str(h5_path), baseline_snapshot)
    roundtrip_snapshot = read_coupling_snapshot(str(h5_path))
    np.testing.assert_allclose(axisymmetric_area_weights(roundtrip_snapshot), area)
    plot_feedback_map(outdir, args, roundtrip_snapshot, r_grid, z_grid)
    plot_feedback_scaling(outdir, rows)
    write_csv(outdir, rows)
    write_summary(outdir, args, rows)
    return rows


def write_csv(outdir: Path, rows: list[dict[str, float | str]]) -> None:
    path = outdir / "axisymmetric_feedback_scan.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_feedback_map(
    outdir: Path,
    args: argparse.Namespace,
    snapshot: CouplingSnapshot,
    current_r_grid: np.ndarray,
    current_z_grid: np.ndarray,
) -> None:
    probe_r_grid, probe_z_grid, probes = probe_grid()
    field = axisymmetric_toroidal_current_field(snapshot, probes)
    feedback = np.linalg.norm(field, axis=1).reshape(probe_r_grid.shape) / args.b_ref
    j_phi = snapshot.toroidal_current_density().reshape(current_r_grid.shape)

    midplane = np.argmin(np.abs(probe_z_grid[:, 0]))
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.8))
    j_img = axes[0].pcolormesh(current_r_grid, current_z_grid, j_phi, shading="auto", cmap="viridis")
    axes[0].set_title(r"manufactured $J_\phi$")
    axes[0].set_xlabel("R [m]")
    axes[0].set_ylabel("Z [m]")
    cbar = fig.colorbar(j_img, ax=axes[0])
    cbar.set_label(r"$J_\phi$ [A/m$^2$]")

    fb_img = axes[1].pcolormesh(
        probe_r_grid,
        probe_z_grid,
        feedback,
        shading="auto",
        cmap="magma",
    )
    axes[1].contour(probe_r_grid, probe_z_grid, feedback, levels=[1.0e-6, 1.0e-5], colors="white", linewidths=0.8)
    axes[1].set_title(r"probe $|\delta B|/|B_{ref}|$")
    axes[1].set_xlabel("R [m]")
    cbar = fig.colorbar(fb_img, ax=axes[1])
    cbar.set_label("feedback ratio")

    axes[2].semilogy(probe_r_grid[midplane, :], feedback[midplane, :], color="C3")
    axes[2].axhline(1.0e-3, color="0.35", ls="--", lw=1.0, label="two-way prompt")
    axes[2].axhline(1.0e-2, color="0.15", ls=":", lw=1.0, label="strong prompt")
    axes[2].set_title("midplane probe line")
    axes[2].set_xlabel("R [m]")
    axes[2].set_ylabel(r"$|\delta B|/|B_{ref}|$")
    axes[2].legend(fontsize=8)
    axes[2].grid(True, which="both", alpha=0.25)

    fig.suptitle(
        f"Axisymmetric current-feedback diagnostic, peak J_phi={args.baseline_j:.1e} A/m^2"
    )
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_feedback_field_map.png", dpi=180)
    plt.close(fig)


def plot_feedback_scaling(outdir: Path, rows: list[dict[str, float | str]]) -> None:
    current = np.array([float(row["absolute_toroidal_current_A"]) for row in rows])
    feedback = np.array([float(row["max_feedback"]) for row in rows])

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(current, feedback, "o-", label="loop-field estimate")
    ax.loglog(current, feedback[0] * current / current[0], "--", color="0.45", label="linear reference")
    ax.axhline(1.0e-3, color="0.35", ls="--", lw=1.0, label="two-way prompt")
    ax.axhline(1.0e-2, color="0.15", ls=":", lw=1.0, label="strong prompt")
    ax.set_xlim(current.min() * 0.75, current.max() * 1.45)
    ax.set_ylim(feedback.min() * 0.7, feedback.max() * 1.6)
    for idx, (xval, yval, row) in enumerate(zip(current, feedback, rows)):
        offset = (-34, -2) if idx == len(rows) - 1 else (4, 4)
        ax.annotate(
            f"{float(row['peak_j_phi_A_m2']):.0e}",
            (xval, yval),
            textcoords="offset points",
            xytext=offset,
            fontsize=7,
        )
    ax.set_xlabel(r"$\int |J_\phi| dA$ [A]")
    ax.set_ylabel(r"max probe $|\delta B|/|B_{ref}|$")
    ax.set_title("Feedback scales linearly with manufactured current")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_feedback_scaling.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, args: argparse.Namespace, rows: list[dict[str, float | str]]) -> None:
    baseline = min(rows, key=lambda row: abs(float(row["peak_j_phi_A_m2"]) - args.baseline_j))
    lines = [
        "# Axisymmetric current-feedback diagnostic",
        "",
        "This is an analytic/reduced-model coupling diagnostic, not a MUG or TokaMaker result.",
        "A manufactured toroidal liquid-metal current patch is converted to circular-loop",
        "current elements and evaluated at plasma-side probe points with the standard",
        "axisymmetric Biot-Savart formula.",
        "",
        "Baseline result:",
        f"- Peak J_phi = {float(baseline['peak_j_phi_A_m2']):.3e} A/m^2.",
        f"- Integrated |J_phi| dA = {float(baseline['absolute_toroidal_current_A']):.3e} A.",
        f"- Max probe |delta B| = {float(baseline['max_delta_B_T']):.3e} T.",
        f"- Max |delta B|/|B_ref| = {float(baseline['max_feedback']):.3e}.",
        f"- Feedback classification = {baseline['classification']}.",
        "",
        "Generated files:",
        "- `axisymmetric_feedback_snapshot.h5`: HDF5 roundtrip snapshot for the baseline current patch.",
        "- `axisymmetric_feedback_scan.csv`: current-amplitude sweep and feedback metrics.",
        "- `axisymmetric_feedback_field_map.png`: current patch, probe feedback map, and midplane line.",
        "- `axisymmetric_feedback_scaling.png`: linearity and threshold plot.",
        "",
        "Red-team interpretation:",
        "- The field estimate includes geometric 1/R loop effects, but no vessel, wall-current,",
        "  iron, plasma-response, or TokaMaker equilibrium solve.",
        "- The sampled currents are manufactured and do not satisfy an electric-potential closure.",
        "- Probes are placed outside the current patch to avoid unresolved self-field singularities.",
        "- Crossing 1e-3 in |delta B|/|B_ref| is treated as a prompt for a coupled sensitivity",
        "  study, not as a prediction of the plasma-equilibrium displacement.",
    ]
    (outdir / "axisymmetric_feedback_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    rows = run_scan(args, outdir)
    print(f"Wrote {len(rows)} axisymmetric feedback diagnostic rows to {outdir}")


if __name__ == "__main__":
    main()
