#!/usr/bin/env python3
"""Audit an LM-MHD CouplingSnapshot HDF5 file.

This is an offline diagnostic for future MUG-to-TokaMaker handoff files. It
computes reduced magnetic-feedback metrics and, when the snapshot points form a
complete structured R-Z grid, conservative axisymmetric current-closure checks.
It does not run MUG or TokaMaker.
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
        axisymmetric_area_weights,
        axisymmetric_boundary_current_integrals,
        axisymmetric_cell_volumes,
        axisymmetric_charge_conservation_residual,
        axisymmetric_net_boundary_current,
        axisymmetric_wall_normal_current_extrema,
        feedback_metrics_from_snapshot,
        read_coupling_snapshot,
        structured_snapshot_current_grid,
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

    axisymmetric_area_weights = coupling.axisymmetric_area_weights
    feedback_metrics_from_snapshot = coupling.feedback_metrics_from_snapshot
    read_coupling_snapshot = coupling.read_coupling_snapshot
    structured_snapshot_current_grid = coupling.structured_snapshot_current_grid
    axisymmetric_boundary_current_integrals = inductionless.axisymmetric_boundary_current_integrals
    axisymmetric_cell_volumes = inductionless.axisymmetric_cell_volumes
    axisymmetric_charge_conservation_residual = inductionless.axisymmetric_charge_conservation_residual
    axisymmetric_net_boundary_current = inductionless.axisymmetric_net_boundary_current
    axisymmetric_wall_normal_current_extrema = inductionless.axisymmetric_wall_normal_current_extrema


DEFAULT_SNAPSHOT = "cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_snapshot.h5"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", default=DEFAULT_SNAPSHOT, help="CouplingSnapshot HDF5 file")
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--prefix", default="coupling_snapshot_audit", help="output filename prefix")
    parser.add_argument("--reference-field", type=float, default=5.0, help="feedback normalization field [T]")
    parser.add_argument("--regularization", type=float, default=0.03, help="loop-field regularization length [m]")
    parser.add_argument(
        "--probe-points",
        default="0.85:0.0,1.10:0.35,1.45:0.0,1.90:-0.35,2.15:0.0",
        help="comma-separated R:Z probe list [m]",
    )
    return parser.parse_args()


def parse_probe_points(raw: str) -> np.ndarray:
    points = []
    for item in raw.split(","):
        if not item.strip():
            continue
        parts = item.split(":")
        if len(parts) != 2:
            raise ValueError("probe points must use R:Z entries")
        points.append((float(parts[0]), float(parts[1])))
    if not points:
        raise ValueError("at least one probe point is required")
    return np.asarray(points, dtype=np.float64)


def audit_snapshot(args: argparse.Namespace) -> tuple[dict[str, float | str], dict[str, np.ndarray] | None]:
    snapshot = read_coupling_snapshot(str(ROOT / args.snapshot if not Path(args.snapshot).is_absolute() else args.snapshot))
    area_weights = axisymmetric_area_weights(snapshot)
    probes = parse_probe_points(args.probe_points)
    metrics = feedback_metrics_from_snapshot(
        snapshot,
        probes,
        reference_field=args.reference_field,
        area_weights=area_weights,
        regularization=args.regularization,
    )
    current_norm = np.linalg.norm(snapshot.current_density, axis=1)
    row: dict[str, float | str] = {
        "snapshot": str(args.snapshot),
        "time": float(snapshot.time),
        "n_points": int(snapshot.points_rz.shape[0]),
        "total_toroidal_current": metrics.total_toroidal_current,
        "absolute_toroidal_current": metrics.absolute_toroidal_current,
        "current_centroid_r": metrics.current_centroid_rz[0],
        "current_centroid_z": metrics.current_centroid_rz[1],
        "max_delta_b": metrics.max_delta_b,
        "rms_delta_b": metrics.rms_delta_b,
        "max_feedback": metrics.max_feedback,
        "rms_feedback": metrics.rms_feedback,
        "classification": metrics.classification,
        "max_abs_j_phi": float(np.max(np.abs(snapshot.current_density[:, 1]))),
        "max_abs_j_poloidal": float(np.max(np.linalg.norm(snapshot.current_density[:, [0, 2]], axis=1))),
        "max_abs_j": float(np.max(current_norm)),
        "structured_grid": "false",
    }
    try:
        r, z, current_grid = structured_snapshot_current_grid(snapshot)
    except ValueError as exc:
        row["structured_grid_reason"] = str(exc)
        return row, None

    residual = axisymmetric_charge_conservation_residual(current_grid, r, z)
    wall_current = axisymmetric_wall_normal_current_extrema(current_grid)
    boundary_current = axisymmetric_boundary_current_integrals(current_grid, r, z)
    residual_region = residual[2:-2, 2:-2] if min(residual.shape) > 4 else residual
    order = np.lexsort((snapshot.points_rz[:, 0], snapshot.points_rz[:, 1]))
    volume_grid = snapshot.cell_volume[order].reshape(z.size, r.size)
    expected_volume = axisymmetric_cell_volumes(r, z)
    volume_scale = max(np.finfo(float).tiny, float(np.max(np.abs(expected_volume))))
    volume_error = float(np.max(np.abs(volume_grid - expected_volume)) / volume_scale)
    row.update(
        {
            "structured_grid": "true",
            "nr": int(r.size),
            "nz": int(z.size),
            "max_cell_volume_rel_error": volume_error,
            "max_abs_axisym_divJ": float(np.max(np.abs(residual_region))),
            "max_wall_normal_current": float(max(wall_current.values())),
            "max_abs_boundary_current": float(max(abs(value) for value in boundary_current.values())),
            "net_boundary_current": float(axisymmetric_net_boundary_current(current_grid, r, z)),
        }
    )
    fields = {
        "r": r,
        "z": z,
        "current": current_grid,
        "residual": residual,
    }
    return row, fields


def write_metrics(outdir: Path, prefix: str, row: dict[str, float | str]) -> None:
    with (outdir / f"{prefix}_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def plot_structured_current(outdir: Path, prefix: str, fields: dict[str, np.ndarray] | None) -> None:
    if fields is None:
        return
    r = fields["r"]
    z = fields["z"]
    current = fields["current"]
    residual = fields["residual"]
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    panels = [
        (current[..., 1], r"$J_\phi$ [A/m$^2$]", "coolwarm"),
        (np.linalg.norm(current[..., [0, 2]], axis=-1), r"$|J_p|$ [A/m$^2$]", "viridis"),
        (np.linalg.norm(current, axis=-1), r"$|J|$ [A/m$^2$]", "magma"),
        (residual, r"axisym. $\nabla\cdot J$", "coolwarm"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    for ax, (values, title, cmap) in zip(axes, panels):
        image = ax.pcolormesh(r_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_title(title)
        fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(outdir / f"{prefix}_current.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, prefix: str, row: dict[str, float | str], plotted: bool) -> None:
    lines = [
        "# CouplingSnapshot audit",
        "",
        "This offline diagnostic reads a `CouplingSnapshot` HDF5 file and",
        "computes reduced current-feedback and current-closure checks. It does",
        "not run MUG or TokaMaker.",
        "",
        "Headline metrics:",
        f"- snapshot = `{row['snapshot']}`",
        f"- points = {int(row['n_points'])}",
        f"- total toroidal current = {float(row['total_toroidal_current']):.3e} A",
        f"- absolute toroidal current = {float(row['absolute_toroidal_current']):.3e} A",
        f"- max |delta B|/|B_ref| = {float(row['max_feedback']):.3e}",
        f"- classification = `{row['classification']}`",
        f"- max |J_phi| = {float(row['max_abs_j_phi']):.3e} A/m^2",
        f"- structured grid = `{row['structured_grid']}`",
    ]
    if row["structured_grid"] == "true":
        lines.extend(
            [
                f"- grid = {int(row['nr'])} x {int(row['nz'])}",
                f"- max cell-volume relative error = {float(row['max_cell_volume_rel_error']):.3e}",
                f"- max interior axisymmetric |div J| = {float(row['max_abs_axisym_divJ']):.3e}",
                f"- max wall-normal current = {float(row['max_wall_normal_current']):.3e} A/m^2",
                f"- max integrated boundary current = {float(row['max_abs_boundary_current']):.3e} A",
                f"- net integrated boundary current = {float(row['net_boundary_current']):.3e} A",
            ]
        )
    else:
        lines.append(f"- structured-grid audit skipped: `{row['structured_grid_reason']}`")
    lines.extend(
        [
            "",
            "Generated files:",
            f"- `{prefix}_metrics.csv`",
        ]
    )
    if plotted:
        lines.append(f"- `{prefix}_current.png`")
    lines.extend(
        [
            "",
            "Red-team interpretation:",
            "- Magnetic feedback is a circular-loop triage estimate.",
            "- Boundary-current checks are only available for complete structured",
            "  R-Z snapshots.",
            "- A clean audit does not validate MUG physics; it verifies this HDF5",
            "  handoff file against reduced postprocessing conventions.",
        ]
    )
    (outdir / f"{prefix}_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    row, fields = audit_snapshot(args)
    write_metrics(outdir, args.prefix, row)
    plot_structured_current(outdir, args.prefix, fields)
    write_summary(outdir, args.prefix, row, plotted=fields is not None)
    print(f"Wrote CouplingSnapshot audit diagnostics to {outdir}")


if __name__ == "__main__":
    main()
