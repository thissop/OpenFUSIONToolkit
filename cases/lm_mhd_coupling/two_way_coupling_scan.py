#!/usr/bin/env python3
"""Regime scan for LM-MHD two-way feedback and unresolved-layer drag triage.

This script is an executable planning aid. It does not run MUG or TokaMaker.
It evaluates standard dimensionless groups over a velocity/field grid and
writes a CSV plus smplotlib-styled figures for red-team review.
"""
from __future__ import annotations

import argparse
import csv
import os
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
        blanket_feedback_ratio_sheet,
        classify_coupling,
        hartmann_number,
        interaction_parameter,
        magnetic_reynolds_number,
        representative_liquid_metals,
        reynolds_number,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    coupling_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = coupling
    spec.loader.exec_module(coupling)
    blanket_feedback_ratio_sheet = coupling.blanket_feedback_ratio_sheet
    classify_coupling = coupling.classify_coupling
    hartmann_number = coupling.hartmann_number
    interaction_parameter = coupling.interaction_parameter
    magnetic_reynolds_number = coupling.magnetic_reynolds_number
    representative_liquid_metals = coupling.representative_liquid_metals
    reynolds_number = coupling.reynolds_number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=float, default=0.05, help="duct half-width / MHD length scale [m]")
    parser.add_argument("--thickness", type=float, default=0.05, help="conducting LM sheet thickness for feedback estimate [m]")
    parser.add_argument("--b-ref", type=float, default=5.0, help="reference field for |delta B|/|B_ref|")
    parser.add_argument(
        "--closure-factor",
        type=float,
        default=0.1,
        help="fraction of sigma U B used in sheet-current upper-bound estimate",
    )
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    return parser.parse_args()


def build_scan(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, dict[str, dict[str, np.ndarray]]]:
    fields = np.geomspace(0.5, 10.0, 80)
    velocities = np.geomspace(1.0e-4, 2.0, 90)
    B_grid, U_grid = np.meshgrid(fields, velocities, indexing="xy")
    data = {}
    for key, props in representative_liquid_metals().items():
        ha = hartmann_number(
            B_grid,
            args.length,
            props.electrical_conductivity,
            props.density,
            props.kinematic_viscosity,
        )
        re = reynolds_number(U_grid, args.length, props.kinematic_viscosity)
        rm = magnetic_reynolds_number(U_grid, args.length, props.electrical_conductivity)
        interaction = interaction_parameter(
            B_grid,
            args.length,
            U_grid,
            props.electrical_conductivity,
            props.density,
        )
        feedback = blanket_feedback_ratio_sheet(
            props,
            U_grid,
            B_grid,
            args.thickness,
            args.b_ref,
            closure_factor=args.closure_factor,
        )
        data[key] = {
            "Ha": ha,
            "Re": re,
            "Rm": rm,
            "N": interaction,
            "feedback": feedback,
        }
    return fields, velocities, data


def write_csv(outdir: Path, fields: np.ndarray, velocities: np.ndarray, data: dict[str, dict[str, np.ndarray]]) -> None:
    path = outdir / "coupling_regime_scan.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "material",
                "B_T",
                "U_m_s",
                "Ha",
                "Re",
                "Rm",
                "N",
                "feedback_sheet_upper",
                "classification",
            ]
        )
        for material, values in data.items():
            for i, U in enumerate(velocities):
                for j, B in enumerate(fields):
                    rm = float(values["Rm"][i, j])
                    feedback = float(values["feedback"][i, j])
                    interaction = float(values["N"][i, j])
                    writer.writerow(
                        [
                            material,
                            f"{B:.8e}",
                            f"{U:.8e}",
                            f"{values['Ha'][i, j]:.8e}",
                            f"{values['Re'][i, j]:.8e}",
                            f"{rm:.8e}",
                            f"{interaction:.8e}",
                            f"{feedback:.8e}",
                            classify_coupling(rm=rm, feedback=feedback, interaction=interaction),
                        ]
                    )


def add_contours(ax: plt.Axes, fields: np.ndarray, velocities: np.ndarray, values: dict[str, np.ndarray]) -> None:
    ax.contour(fields, velocities, values["N"], levels=[1.0, 10.0, 100.0], colors="white", linewidths=0.8)
    ax.contour(fields, velocities, values["Rm"], levels=[0.01, 0.1], colors="k", linestyles="--", linewidths=0.8)
    ax.contour(fields, velocities, values["feedback"], levels=[1.0e-3, 1.0e-2], colors="0.2", linewidths=0.9)


def plot_heatmaps(outdir: Path, fields: np.ndarray, velocities: np.ndarray, data: dict[str, dict[str, np.ndarray]]) -> None:
    for material, values in data.items():
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)

        n_img = axes[0].pcolormesh(fields, velocities, values["N"], norm="log", shading="auto", cmap="viridis")
        axes[0].set_title(f"{material}: interaction parameter N")
        axes[0].set_xlabel("B [T]")
        axes[0].set_ylabel("U [m/s]")
        cbar = fig.colorbar(n_img, ax=axes[0])
        cbar.set_label("N")
        add_contours(axes[0], fields, velocities, values)

        fb_img = axes[1].pcolormesh(
            fields,
            velocities,
            values["feedback"],
            norm="log",
            shading="auto",
            cmap="magma",
        )
        axes[1].set_title(f"{material}: |delta B|/|B_ref| sheet estimate")
        axes[1].set_xlabel("B [T]")
        cbar = fig.colorbar(fb_img, ax=axes[1])
        cbar.set_label("upper-bound feedback ratio")
        add_contours(axes[1], fields, velocities, values)

        for ax in axes:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.25)

        fig.suptitle("Contours: N={1,10,100}, Rm={0.01,0.1} dashed, feedback={1e-3,1e-2}")
        fig.tight_layout()
        fig.savefig(outdir / f"{material.lower()}_regime_heatmap.png", dpi=180)
        plt.close(fig)


def plot_reduced_lines(outdir: Path, fields: np.ndarray, velocities: np.ndarray, data: dict[str, dict[str, np.ndarray]]) -> None:
    b_target = 5.0
    b_idx = int(np.argmin(np.abs(fields - b_target)))
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax2 = ax.twinx()
    for material, values in data.items():
        ax.loglog(velocities, values["Rm"][:, b_idx], label=f"{material} Rm")
        ax2.loglog(velocities, values["feedback"][:, b_idx], "--", label=f"{material} feedback")
    ax.axhline(0.1, color="0.3", ls=":", lw=1.0)
    ax2.axhline(1.0e-3, color="0.5", ls="-.", lw=1.0)
    ax.set_xlabel("U [m/s]")
    ax.set_ylabel("Rm")
    ax2.set_ylabel("upper-bound |delta B|/|B_ref|")
    ax.set_title(f"Two-way triage at B ~= {fields[b_idx]:.2f} T")
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc="upper left", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "two_way_threshold_lines.png", dpi=180)
    plt.close(fig)


def write_summary(
    outdir: Path,
    args: argparse.Namespace,
    fields: np.ndarray,
    velocities: np.ndarray,
    data: dict[str, dict[str, np.ndarray]],
) -> None:
    lines = [
        "# LM-MHD coupling regime scan",
        "",
        "This is an executable triage artifact, not a MUG/TokaMaker result. It scans representative",
        "liquid-metal properties and reports dimensionless groups plus an intentionally conservative",
        "sheet-current feedback estimate.",
        "",
        "Assumptions:",
        f"- Length scale L = {args.length:g} m.",
        f"- Conducting sheet thickness for feedback estimate = {args.thickness:g} m.",
        f"- Reference field for feedback ratio = {args.b_ref:g} T.",
        f"- Closure factor multiplying sigma U B = {args.closure_factor:g}. A value below one",
        "  acknowledges that electric-potential closure can reduce or redirect the current.",
        "- Classification thresholds: Rm >= 0.1 => question inductionless; feedback >= 1e-3 =>",
        "  two-way-candidate; otherwise N >= 1 => drag-important.",
        "",
        "| material | Ha range | Re range | Rm range | N range | feedback range | dominant labels |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for material, values in data.items():
        labels = {}
        for i in range(values["Rm"].shape[0]):
            for j in range(values["Rm"].shape[1]):
                label = classify_coupling(
                    rm=float(values["Rm"][i, j]),
                    feedback=float(values["feedback"][i, j]),
                    interaction=float(values["N"][i, j]),
                )
                labels[label] = labels.get(label, 0) + 1
        label_text = ", ".join(f"{key}:{value}" for key, value in sorted(labels.items()))
        lines.append(
            f"| {material} | {values['Ha'].min():.2e}-{values['Ha'].max():.2e} | "
            f"{values['Re'].min():.2e}-{values['Re'].max():.2e} | "
            f"{values['Rm'].min():.2e}-{values['Rm'].max():.2e} | "
            f"{values['N'].min():.2e}-{values['N'].max():.2e} | "
            f"{values['feedback'].min():.2e}-{values['feedback'].max():.2e} | {label_text} |"
        )
    lines += [
        "",
        "Red-team interpretation:",
        "- Large N with small Rm means one-way imposed-field LM-MHD can still have very strong",
        "  Lorentz braking and unresolved Hartmann-layer pressure loss. This is the regime where",
        "  effective drag/SM82-style closures are most valuable for full-device meshes.",
        "- A feedback ratio above 1e-3 is a prompt to run a coupled field sensitivity study, not",
        "  proof that the plasma equilibrium will change by that amount.",
        "- The sheet-current estimate is an upper-bound scale. A real duct can have smaller or",
        "  differently directed currents after solving electric potential and wall closure.",
        "",
        "Generated files:",
        "- `coupling_regime_scan.csv`",
        "- `pbli_regime_heatmap.png`",
        "- `li_regime_heatmap.png`",
        "- `two_way_threshold_lines.png`",
    ]
    (outdir / "summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    fields, velocities, data = build_scan(args)
    write_csv(outdir, fields, velocities, data)
    plot_heatmaps(outdir, fields, velocities, data)
    plot_reduced_lines(outdir, fields, velocities, data)
    write_summary(outdir, args, fields, velocities, data)
    print(f"Wrote LM-MHD coupling scan artifacts to {outdir}")


if __name__ == "__main__":
    main()
