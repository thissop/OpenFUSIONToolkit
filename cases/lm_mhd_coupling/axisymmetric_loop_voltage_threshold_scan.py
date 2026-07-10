#!/usr/bin/env python3
"""Loop-voltage/conductivity threshold scan for axisymmetric feedback.

This wraps the manufactured loop-voltage feedback diagnostic over a small
parameter grid to estimate when max |delta B|/|B_ref| crosses the planning
thresholds for two-way coupling. It is a reduced diagnostic, not a coupled
equilibrium solve.
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
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from axisymmetric_loop_voltage_feedback import metric_row  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--nr", type=int, default=41, help="radial grid points")
    parser.add_argument("--nz", type=int, default=33, help="vertical grid points")
    parser.add_argument("--b-phi", type=float, default=5.0, help="dominant toroidal magnetic field [T]")
    parser.add_argument("--reference-field", type=float, default=5.0, help="feedback normalization field [T]")
    parser.add_argument("--sigma-values", default="2e5,5e5,8e5,1.2e6,2e6", help="conductivity scan [S/m]")
    parser.add_argument("--loop-voltage-values", default="0,0.05,0.1,0.2,0.5,1.0", help="loop-voltage scan [V]")
    return parser.parse_args()


def scan(args: argparse.Namespace) -> tuple[list[dict[str, float | str]], np.ndarray, np.ndarray, np.ndarray]:
    sigma_values = np.array([float(item) for item in args.sigma_values.split(",") if item.strip()])
    voltages = np.array([float(item) for item in args.loop_voltage_values.split(",") if item.strip()])
    feedback = np.zeros((sigma_values.size, voltages.size), dtype=np.float64)
    rows: list[dict[str, float | str]] = []
    for i, sigma0 in enumerate(sigma_values):
        local_args = argparse.Namespace(
            outdir=args.outdir,
            nr=args.nr,
            nz=args.nz,
            sigma0=float(sigma0),
            b_phi=args.b_phi,
            reference_field=args.reference_field,
        )
        for j, voltage in enumerate(voltages):
            row = metric_row(local_args, float(voltage))
            row["sigma0"] = float(sigma0)
            rows.append(row)
            feedback[i, j] = float(row["max_feedback"])
    return rows, sigma_values, voltages, feedback


def threshold_voltage(voltages: np.ndarray, feedback: np.ndarray, threshold: float) -> float:
    if np.max(feedback) < threshold:
        return float("nan")
    if np.min(feedback) >= threshold:
        return float(voltages[0])
    order = np.argsort(voltages)
    return float(np.interp(threshold, feedback[order], voltages[order]))


def threshold_table(sigma_values: np.ndarray, voltages: np.ndarray, feedback: np.ndarray) -> list[dict[str, float]]:
    rows = []
    for sigma0, values in zip(sigma_values, feedback):
        rows.append(
            {
                "sigma0": float(sigma0),
                "vloop_two_way_threshold": threshold_voltage(voltages, values, 1.0e-3),
                "vloop_strong_threshold": threshold_voltage(voltages, values, 1.0e-2),
                "max_feedback_at_max_voltage": float(values[-1]),
            }
        )
    return rows


def write_csvs(outdir: Path, rows: list[dict[str, float | str]], thresholds: list[dict[str, float]]) -> None:
    with (outdir / "axisymmetric_loop_voltage_threshold_scan.csv").open("w", newline="") as handle:
        keys = ["sigma0"] + [key for key in rows[0] if key != "sigma0"]
        writer = csv.DictWriter(handle, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    with (outdir / "axisymmetric_loop_voltage_thresholds.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(thresholds[0]), lineterminator="\n")
        writer.writeheader()
        for row in thresholds:
            writer.writerow(row)


def plot_heatmap(outdir: Path, sigma_values: np.ndarray, voltages: np.ndarray, feedback: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(7.3, 4.8))
    image = ax.pcolormesh(voltages, sigma_values, feedback, shading="auto", cmap="viridis")
    contours = ax.contour(voltages, sigma_values, feedback, levels=[1.0e-3, 1.0e-2], colors=["white", "red"])
    ax.clabel(contours, fmt={1.0e-3: "1e-3", 1.0e-2: "1e-2"}, fontsize=8)
    ax.set_xlabel(r"$V_{loop}$ [V]")
    ax.set_ylabel(r"$\sigma_0$ [S/m]")
    ax.set_title(r"max $|\delta B|/|B_{ref}|$")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_loop_voltage_threshold_heatmap.png", dpi=180)
    plt.close(fig)


def plot_thresholds(outdir: Path, thresholds: list[dict[str, float]]) -> None:
    sigma = np.array([row["sigma0"] for row in thresholds])
    two_way = np.array([row["vloop_two_way_threshold"] for row in thresholds])
    strong = np.array([row["vloop_strong_threshold"] for row in thresholds])
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(sigma, two_way, "o-", label=r"$|\delta B|/B=10^{-3}$")
    ax.plot(sigma, strong, "s--", label=r"$|\delta B|/B=10^{-2}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sigma_0$ [S/m]")
    ax.set_ylabel(r"threshold $V_{loop}$ [V]")
    ax.set_title("Loop voltage threshold estimate")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "axisymmetric_loop_voltage_thresholds.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, thresholds: list[dict[str, float]]) -> None:
    pblike = min(thresholds, key=lambda row: abs(row["sigma0"] - 8.0e5))
    pblike_strong = "not reached" if np.isnan(pblike["vloop_strong_threshold"]) else f"{pblike['vloop_strong_threshold']:.3e} V"
    lines = [
        "# Axisymmetric loop-voltage/conductivity threshold scan",
        "",
        "This reduced manufactured scan estimates when loop-voltage-driven",
        "axisymmetric blanket currents cross magnetic-feedback planning",
        "thresholds. It does not run MUG or TokaMaker.",
        "",
        "PbLi-like row nearest `sigma0 = 8e5 S/m`:",
        f"- two-way threshold (`1e-3`) V_loop = {pblike['vloop_two_way_threshold']:.3e} V.",
        f"- strong threshold (`1e-2`) V_loop = {pblike_strong}.",
        f"- max feedback at max scanned voltage = {pblike['max_feedback_at_max_voltage']:.3e}.",
        "",
        "Threshold table:",
        "",
        "| sigma0 [S/m] | V_loop at 1e-3 [V] | V_loop at 1e-2 [V] | max feedback |",
        "|---:|---:|---:|---:|",
    ]
    for row in thresholds:
        strong = "not reached" if np.isnan(row["vloop_strong_threshold"]) else f"`{row['vloop_strong_threshold']:.3e}`"
        two_way = "not reached" if np.isnan(row["vloop_two_way_threshold"]) else f"`{row['vloop_two_way_threshold']:.3e}`"
        lines.append(
            f"| `{row['sigma0']:.3e}` | {two_way} | {strong} | `{row['max_feedback_at_max_voltage']:.3e}` |"
        )
    lines.extend(
        [
            "",
            "Generated files:",
            "- `axisymmetric_loop_voltage_threshold_scan.csv`",
            "- `axisymmetric_loop_voltage_thresholds.csv`",
            "- `axisymmetric_loop_voltage_threshold_heatmap.png`",
            "- `axisymmetric_loop_voltage_thresholds.png`",
            "",
            "Red-team interpretation:",
            "- Thresholds are interpolation over a coarse manufactured scan.",
            "- The current closure, geometry, and field feedback are reduced models.",
            "- Use these numbers only to prioritize coupled sensitivity studies, not",
            "  as blanket design predictions.",
        ]
    )
    (outdir / "axisymmetric_loop_voltage_threshold_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    rows, sigma_values, voltages, feedback = scan(args)
    thresholds = threshold_table(sigma_values, voltages, feedback)
    write_csvs(outdir, rows, thresholds)
    plot_heatmap(outdir, sigma_values, voltages, feedback)
    plot_thresholds(outdir, thresholds)
    write_summary(outdir, thresholds)
    print(f"Wrote axisymmetric loop-voltage threshold diagnostics to {outdir}")


if __name__ == "__main__":
    main()
