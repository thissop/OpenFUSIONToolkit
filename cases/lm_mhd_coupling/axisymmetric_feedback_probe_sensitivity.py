#!/usr/bin/env python3
"""Probe-location sensitivity for reduced axisymmetric feedback metrics.

This red-team diagnostic asks whether a CouplingSnapshot's two-way feedback
classification is robust to probe placement and circular-loop regularization.
It is a reduced postprocessing scan, not a TokaMaker response calculation.
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
        feedback_metrics_from_snapshot,
        read_coupling_snapshot,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    coupling_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    coupling_spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(coupling_spec)
    sys.modules[coupling_spec.name] = coupling
    coupling_spec.loader.exec_module(coupling)
    axisymmetric_area_weights = coupling.axisymmetric_area_weights
    feedback_metrics_from_snapshot = coupling.feedback_metrics_from_snapshot
    read_coupling_snapshot = coupling.read_coupling_snapshot


DEFAULT_SNAPSHOT = "cases/lm_mhd_coupling/results/axisymmetric_loop_voltage_feedback_snapshot.h5"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", default=DEFAULT_SNAPSHOT, help="CouplingSnapshot HDF5 file")
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--reference-field", type=float, default=5.0, help="feedback normalization field [T]")
    parser.add_argument("--offsets", default="0.05,0.1,0.2,0.35,0.5", help="probe offsets from current annulus [m]")
    parser.add_argument(
        "--regularizations",
        default="0.01,0.02,0.03,0.05,0.1",
        help="loop-field regularization lengths [m]",
    )
    parser.add_argument("--z-probes", default="-0.4,-0.2,0.0,0.2,0.4", help="probe Z locations [m]")
    return parser.parse_args()


def snapshot_path(raw: str) -> Path:
    path = Path(raw)
    return path if path.is_absolute() else ROOT / path


def probe_line(side: str, radius: float, z_values: np.ndarray) -> np.ndarray:
    if side not in {"inboard", "outboard"}:
        raise ValueError("side must be 'inboard' or 'outboard'")
    return np.column_stack((np.full_like(z_values, radius, dtype=np.float64), z_values))


def scan(args: argparse.Namespace) -> tuple[list[dict[str, float | str]], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    snapshot = read_coupling_snapshot(str(snapshot_path(args.snapshot)))
    area_weights = axisymmetric_area_weights(snapshot)
    offsets = np.array([float(item) for item in args.offsets.split(",") if item.strip()], dtype=np.float64)
    regularizations = np.array([float(item) for item in args.regularizations.split(",") if item.strip()], dtype=np.float64)
    z_values = np.array([float(item) for item in args.z_probes.split(",") if item.strip()], dtype=np.float64)
    if np.any(offsets <= 0.0):
        raise ValueError("offsets must be positive")
    if np.any(regularizations < 0.0):
        raise ValueError("regularizations must be nonnegative")
    r_min = float(np.min(snapshot.points_rz[:, 0]))
    r_max = float(np.max(snapshot.points_rz[:, 0]))
    feedback = np.zeros((2, regularizations.size, offsets.size), dtype=np.float64)
    rows: list[dict[str, float | str]] = []
    for side_index, side in enumerate(("inboard", "outboard")):
        for reg_index, regularization in enumerate(regularizations):
            for offset_index, offset in enumerate(offsets):
                radius = r_min - offset if side == "inboard" else r_max + offset
                probes = probe_line(side, radius, z_values)
                metrics = feedback_metrics_from_snapshot(
                    snapshot,
                    probes,
                    reference_field=args.reference_field,
                    area_weights=area_weights,
                    regularization=float(regularization),
                )
                feedback[side_index, reg_index, offset_index] = metrics.max_feedback
                rows.append(
                    {
                        "side": side,
                        "offset": float(offset),
                        "probe_radius": float(radius),
                        "regularization": float(regularization),
                        "max_feedback": metrics.max_feedback,
                        "rms_feedback": metrics.rms_feedback,
                        "max_delta_b": metrics.max_delta_b,
                        "classification": metrics.classification,
                    }
                )
    return rows, offsets, regularizations, z_values, feedback


def write_csv(outdir: Path, rows: list[dict[str, float | str]]) -> None:
    with (outdir / "axisymmetric_feedback_probe_sensitivity.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_heatmaps(outdir: Path, offsets: np.ndarray, regularizations: np.ndarray, feedback: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharey=True)
    for ax, side, values in zip(axes, ("inboard", "outboard"), feedback):
        image = ax.pcolormesh(offsets, regularizations, values, shading="auto", cmap="viridis")
        contours = ax.contour(offsets, regularizations, values, levels=[1.0e-3, 1.0e-2], colors=["white", "red"])
        ax.clabel(contours, fmt={1.0e-3: "1e-3", 1.0e-2: "1e-2"}, fontsize=8)
        ax.set_xlabel("probe offset [m]")
        ax.set_title(side)
        ax.grid(True, alpha=0.18)
        fig.colorbar(image, ax=ax)
    axes[0].set_ylabel("regularization [m]")
    fig.suptitle(r"Probe sensitivity of max $|\delta B|/|B_{ref}|$")
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_feedback_probe_sensitivity_heatmap.png", dpi=180)
    plt.close(fig)


def plot_offset_curves(outdir: Path, offsets: np.ndarray, regularizations: np.ndarray, feedback: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.1), sharey=True)
    for ax, side, values in zip(axes, ("inboard", "outboard"), feedback):
        for regularization, curve in zip(regularizations, values):
            ax.semilogy(offsets, curve, "o-", label=f"eps={regularization:g} m")
        ax.axhline(1.0e-3, color="0.35", ls="--", label="1e-3")
        ax.axhline(1.0e-2, color="0.20", ls=":", label="1e-2")
        ax.set_xlabel("probe offset [m]")
        ax.set_title(side)
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(r"max $|\delta B|/|B_{ref}|$")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_feedback_probe_sensitivity_curves.png", dpi=180)
    plt.close(fig)


def summarize_rows(rows: list[dict[str, float | str]]) -> dict[str, float | str]:
    strongest = max(rows, key=lambda row: float(row["max_feedback"]))
    weakest = min(rows, key=lambda row: float(row["max_feedback"]))
    default_like = min(
        rows,
        key=lambda row: abs(float(row["offset"]) - 0.1) + abs(float(row["regularization"]) - 0.03),
    )
    two_way_count = sum(float(row["max_feedback"]) >= 1.0e-3 for row in rows)
    strong_count = sum(float(row["max_feedback"]) >= 1.0e-2 for row in rows)
    return {
        "strongest": strongest,
        "weakest": weakest,
        "default_like": default_like,
        "two_way_fraction": two_way_count / len(rows),
        "strong_fraction": strong_count / len(rows),
    }


def write_summary(outdir: Path, rows: list[dict[str, float | str]]) -> None:
    summary = summarize_rows(rows)
    strongest = summary["strongest"]
    weakest = summary["weakest"]
    default_like = summary["default_like"]
    lines = [
        "# Axisymmetric feedback probe-sensitivity scan",
        "",
        "This red-team diagnostic varies probe placement and circular-loop",
        "regularization for a fixed `CouplingSnapshot`. It does not run MUG or",
        "TokaMaker.",
        "",
        "Default-like row nearest offset `0.1 m` and regularization `0.03 m`:",
        f"- side = `{default_like['side']}`.",
        f"- max |delta B|/|B_ref| = {float(default_like['max_feedback']):.3e}.",
        f"- classification = `{default_like['classification']}`.",
        "",
        "Extremes across the scan:",
        f"- strongest feedback = {float(strongest['max_feedback']):.3e} "
        f"(`{strongest['side']}`, offset {float(strongest['offset']):.3g} m, "
        f"eps {float(strongest['regularization']):.3g} m).",
        f"- weakest feedback = {float(weakest['max_feedback']):.3e} "
        f"(`{weakest['side']}`, offset {float(weakest['offset']):.3g} m, "
        f"eps {float(weakest['regularization']):.3g} m).",
        f"- fraction of rows at or above `1e-3` = {float(summary['two_way_fraction']):.3f}.",
        f"- fraction of rows at or above `1e-2` = {float(summary['strong_fraction']):.3f}.",
        "",
        "Generated files:",
        "- `axisymmetric_feedback_probe_sensitivity.csv`",
        "- `axisymmetric_feedback_probe_sensitivity_heatmap.png`",
        "- `axisymmetric_feedback_probe_sensitivity_curves.png`",
        "",
        "Red-team interpretation:",
        "- This scan tests sensitivity of the reduced circular-loop metric, not",
        "  plasma equilibrium sensitivity.",
        "- Rows close to the current annulus are intentionally near-field-biased.",
        "- If only near-field rows cross `1e-3`, that weakens the two-way case.",
        "  If many offset/regularization rows cross `1e-3`, that strengthens the",
        "  priority for coupled sensitivity studies.",
    ]
    (outdir / "axisymmetric_feedback_probe_sensitivity_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    rows, offsets, regularizations, _, feedback = scan(args)
    write_csv(outdir, rows)
    plot_heatmaps(outdir, offsets, regularizations, feedback)
    plot_offset_curves(outdir, offsets, regularizations, feedback)
    write_summary(outdir, rows)
    print(f"Wrote axisymmetric feedback probe-sensitivity diagnostics to {outdir}")


if __name__ == "__main__":
    main()
