#!/usr/bin/env python3
"""Material B-U regime maps for LM-MHD priority triage.

This reduced scan maps Hartmann number, interaction parameter, magnetic
Reynolds number, and sheet-current feedback estimates over magnetic-field and
velocity ranges for representative liquid metals. It is a planning diagnostic,
not a substitute for resolved MUG, inductionless solves, or TokaMaker response.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
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
    coupling_spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(coupling_spec)
    sys.modules[coupling_spec.name] = coupling
    coupling_spec.loader.exec_module(coupling)
    blanket_feedback_ratio_sheet = coupling.blanket_feedback_ratio_sheet
    classify_coupling = coupling.classify_coupling
    hartmann_number = coupling.hartmann_number
    interaction_parameter = coupling.interaction_parameter
    magnetic_reynolds_number = coupling.magnetic_reynolds_number
    representative_liquid_metals = coupling.representative_liquid_metals
    reynolds_number = coupling.reynolds_number


CLASS_VALUE = {
    "low": 0,
    "drag-important": 1,
    "two-way-candidate": 2,
    "full-induction-candidate": 3,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--length", type=float, default=0.05, help="representative channel/blanket length [m]")
    parser.add_argument("--thickness", type=float, default=0.02, help="conducting sheet thickness for feedback [m]")
    parser.add_argument("--reference-field", type=float, default=5.0, help="feedback normalization field [T]")
    parser.add_argument("--closure-factor", type=float, default=0.25, help="sheet-current closure factor")
    parser.add_argument("--velocity-min", type=float, default=0.01, help="min velocity [m/s]")
    parser.add_argument("--velocity-max", type=float, default=1.0, help="max velocity [m/s]")
    parser.add_argument("--field-min", type=float, default=0.5, help="min magnetic field [T]")
    parser.add_argument("--field-max", type=float, default=10.0, help="max magnetic field [T]")
    parser.add_argument("--n", type=int, default=80, help="points per velocity/field axis")
    return parser.parse_args()


def compute_maps(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, dict[str, dict[str, np.ndarray]]]:
    velocities = np.geomspace(args.velocity_min, args.velocity_max, args.n)
    fields = np.geomspace(args.field_min, args.field_max, args.n)
    b_grid, u_grid = np.meshgrid(fields, velocities, indexing="xy")
    materials = representative_liquid_metals()
    maps: dict[str, dict[str, np.ndarray]] = {}
    for name, props in materials.items():
        re = reynolds_number(u_grid, args.length, props.kinematic_viscosity)
        ha = hartmann_number(
            b_grid,
            args.length,
            props.electrical_conductivity,
            props.density,
            props.kinematic_viscosity,
        )
        interaction = interaction_parameter(
            b_grid,
            args.length,
            u_grid,
            props.electrical_conductivity,
            props.density,
        )
        rm = magnetic_reynolds_number(u_grid, args.length, props.electrical_conductivity)
        feedback = blanket_feedback_ratio_sheet(
            props,
            velocity=u_grid,
            magnetic_field=b_grid,
            conducting_thickness=args.thickness,
            reference_field=args.reference_field,
            closure_factor=args.closure_factor,
        )
        class_grid = np.empty(u_grid.shape, dtype=object)
        class_value = np.zeros(u_grid.shape, dtype=np.float64)
        for index in np.ndindex(u_grid.shape):
            label = classify_coupling(
                rm=float(rm[index]),
                feedback=float(feedback[index]),
                interaction=float(interaction[index]),
            )
            class_grid[index] = label
            class_value[index] = CLASS_VALUE[label]
        maps[name] = {
            "re": re,
            "ha": ha,
            "interaction": interaction,
            "rm": rm,
            "feedback": feedback,
            "class_grid": class_grid,
            "class_value": class_value,
        }
    return velocities, fields, maps


def write_scan_csv(outdir: Path, velocities: np.ndarray, fields: np.ndarray, maps: dict[str, dict[str, np.ndarray]]) -> None:
    with (outdir / "lm_mhd_regime_priority_map.csv").open("w", newline="") as handle:
        keys = ["material", "velocity", "magnetic_field", "re", "ha", "interaction", "rm", "feedback", "classification"]
        writer = csv.DictWriter(handle, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        for material, values in maps.items():
            for i, velocity in enumerate(velocities):
                for j, magnetic_field in enumerate(fields):
                    writer.writerow(
                        {
                            "material": material,
                            "velocity": float(velocity),
                            "magnetic_field": float(magnetic_field),
                            "re": float(values["re"][i, j]),
                            "ha": float(values["ha"][i, j]),
                            "interaction": float(values["interaction"][i, j]),
                            "rm": float(values["rm"][i, j]),
                            "feedback": float(values["feedback"][i, j]),
                            "classification": values["class_grid"][i, j],
                        }
                    )


def plot_metric_maps(outdir: Path, velocities: np.ndarray, fields: np.ndarray, maps: dict[str, dict[str, np.ndarray]]) -> None:
    materials = list(maps)
    fig, axes = plt.subplots(len(materials), 3, figsize=(13.6, 4.4 * len(materials)), sharex=True, sharey=True)
    if len(materials) == 1:
        axes = axes[None, :]
    for row_index, material in enumerate(materials):
        values = maps[material]
        panels = [
            ("interaction", r"$\log_{10} N$", [1.0]),
            ("rm", r"$\log_{10} R_m$", [0.1]),
            ("feedback", r"$\log_{10} |\delta B|/B$", [1.0e-3, 1.0e-2]),
        ]
        for ax, (key, title, levels) in zip(axes[row_index], panels):
            data = values[key]
            image = ax.pcolormesh(fields, velocities, np.log10(data), shading="auto", cmap="viridis")
            contours = ax.contour(fields, velocities, data, levels=levels, colors="white")
            labels = {level: f"{level:g}" for level in levels}
            ax.clabel(contours, fmt=labels, fontsize=8)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("B [T]")
            ax.set_ylabel("U [m/s]")
            ax.set_title(f"{material}: {title}")
            fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(outdir / "lm_mhd_regime_priority_metrics.png", dpi=180)
    plt.close(fig)


def plot_classification_maps(
    outdir: Path,
    velocities: np.ndarray,
    fields: np.ndarray,
    maps: dict[str, dict[str, np.ndarray]],
) -> None:
    materials = list(maps)
    fig, axes = plt.subplots(1, len(materials), figsize=(6.2 * len(materials), 4.8), sharex=True, sharey=True)
    if len(materials) == 1:
        axes = [axes]
    cmap = plt.get_cmap("viridis", len(CLASS_VALUE))
    for ax, material in zip(axes, materials):
        image = ax.pcolormesh(fields, velocities, maps[material]["class_value"], shading="auto", cmap=cmap, vmin=0, vmax=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("B [T]")
        ax.set_ylabel("U [m/s]")
        ax.set_title(material)
        cbar = fig.colorbar(image, ax=ax, ticks=list(CLASS_VALUE.values()))
        cbar.ax.set_yticklabels(list(CLASS_VALUE))
    fig.suptitle("Reduced LM-MHD priority classification")
    fig.tight_layout()
    fig.savefig(outdir / "lm_mhd_regime_priority_classification.png", dpi=180)
    plt.close(fig)


def summarize(maps: dict[str, dict[str, np.ndarray]], velocities: np.ndarray, fields: np.ndarray) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    target_velocity = 0.2
    target_field = 5.0
    i = int(np.argmin(np.abs(velocities - target_velocity)))
    j = int(np.argmin(np.abs(fields - target_field)))
    for material, values in maps.items():
        counts = Counter(values["class_grid"].ravel())
        total = values["class_grid"].size
        row: dict[str, float | str] = {
            "material": material,
            "representative_velocity": float(velocities[i]),
            "representative_field": float(fields[j]),
            "representative_ha": float(values["ha"][i, j]),
            "representative_re": float(values["re"][i, j]),
            "representative_interaction": float(values["interaction"][i, j]),
            "representative_rm": float(values["rm"][i, j]),
            "representative_feedback": float(values["feedback"][i, j]),
            "representative_classification": values["class_grid"][i, j],
        }
        for label in CLASS_VALUE:
            row[f"fraction_{label}"] = counts[label] / total
        rows.append(row)
    return rows


def write_summary(outdir: Path, args: argparse.Namespace, summary_rows: list[dict[str, float | str]]) -> None:
    with (outdir / "lm_mhd_regime_priority_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)
    lines = [
        "# LM-MHD B-U regime priority map",
        "",
        "This reduced map scans representative liquid-metal properties over",
        "magnetic-field and velocity ranges. It is a triage calculation, not a",
        "resolved MUG or TokaMaker result.",
        "",
        "Assumptions:",
        f"- length = {args.length:.3e} m",
        f"- conducting thickness = {args.thickness:.3e} m",
        f"- feedback reference field = {args.reference_field:.3e} T",
        f"- sheet-current closure factor = {args.closure_factor:.3e}",
        "",
        "Representative row nearest `U = 0.2 m/s`, `B = 5 T`:",
        "",
        "| material | Ha | Re | N | Rm | feedback | class |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['material']} | `{row['representative_ha']:.3e}` | "
            f"`{row['representative_re']:.3e}` | `{row['representative_interaction']:.3e}` | "
            f"`{row['representative_rm']:.3e}` | `{row['representative_feedback']:.3e}` | "
            f"`{row['representative_classification']}` |"
        )
    lines.extend(
        [
            "",
            "Classification fractions over the scanned B-U grid:",
            "",
            "| material | low | drag-important | two-way-candidate | full-induction-candidate |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in summary_rows:
        lines.append(
            f"| {row['material']} | `{row['fraction_low']:.3f}` | "
            f"`{row['fraction_drag-important']:.3f}` | `{row['fraction_two-way-candidate']:.3f}` | "
            f"`{row['fraction_full-induction-candidate']:.3f}` |"
        )
    lines.extend(
        [
            "",
            "Generated files:",
            "- `lm_mhd_regime_priority_map.csv`",
            "- `lm_mhd_regime_priority_summary.csv`",
            "- `lm_mhd_regime_priority_metrics.png`",
            "- `lm_mhd_regime_priority_classification.png`",
            "",
            "Red-team interpretation:",
            "- `N` flags Lorentz-force/drag importance; it says nothing by itself",
            "  about magnetic feedback into the plasma.",
            "- `Rm` flags when inductionless assumptions need scrutiny.",
            "- The feedback panel is a sheet-current scale with a closure factor,",
            "  not a solved current-closure or equilibrium perturbation.",
            "- These maps prioritize what to simulate next; they do not validate",
            "  MUG or claim blanket design performance.",
        ]
    )
    (outdir / "lm_mhd_regime_priority_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    velocities, fields, maps = compute_maps(args)
    write_scan_csv(outdir, velocities, fields, maps)
    plot_metric_maps(outdir, velocities, fields, maps)
    plot_classification_maps(outdir, velocities, fields, maps)
    write_summary(outdir, args, summarize(maps, velocities, fields))
    print(f"Wrote LM-MHD regime priority diagnostics to {outdir}")


if __name__ == "__main__":
    main()
