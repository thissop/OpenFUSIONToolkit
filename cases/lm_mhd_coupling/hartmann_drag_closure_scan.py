#!/usr/bin/env python3
"""Effective Hartmann-drag closure diagnostics for underresolved LM meshes.

This script is an analytic/reduced-model planning artifact. It does not run MUG
or TokaMaker. It saves CSV/PNG outputs that make the resolved-vs-closure trade
visible for full-device liquid-metal blanket meshes.
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
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "python"))

try:
    from OpenFUSIONToolkit.LM_MHD import (  # noqa: E402
        CouplingSnapshot,
        cells_per_hartmann_layer,
        classify_coupling,
        classify_hartmann_resolution,
        feedback_ratio,
        hartmann_channel_effective_drag,
        hartmann_layer_thickness,
        read_coupling_snapshot,
        sheet_delta_b,
        write_coupling_snapshot,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    def load_module(name: str, relpath: str):
        path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / relpath
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module

    closures = load_module("lm_mhd_closures", "closures.py")
    coupling = load_module("lm_mhd_coupling", "coupling.py")
    CouplingSnapshot = coupling.CouplingSnapshot
    cells_per_hartmann_layer = closures.cells_per_hartmann_layer
    classify_coupling = coupling.classify_coupling
    classify_hartmann_resolution = closures.classify_hartmann_resolution
    feedback_ratio = coupling.feedback_ratio
    hartmann_channel_effective_drag = closures.hartmann_channel_effective_drag
    hartmann_layer_thickness = closures.hartmann_layer_thickness
    read_coupling_snapshot = coupling.read_coupling_snapshot
    sheet_delta_b = coupling.sheet_delta_b
    write_coupling_snapshot = coupling.write_coupling_snapshot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--half-width", type=float, default=0.05, help="channel half-width a [m]")
    parser.add_argument("--nu", type=float, default=1.8e-7, help="kinematic viscosity [m^2/s]")
    parser.add_argument("--target-cells", type=float, default=8.0, help="cells required across a Hartmann layer")
    parser.add_argument("--thickness", type=float, default=0.05, help="LM sheet thickness for current-feedback diagnostic [m]")
    parser.add_argument("--b-ref", type=float, default=5.0, help="reference field for feedback diagnostic [T]")
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    return parser.parse_args()


def build_hartmann_scan(args: argparse.Namespace) -> dict[str, np.ndarray]:
    ha = np.geomspace(1.0, 500.0, 120)
    cells_half = np.unique(np.rint(np.geomspace(4.0, 2000.0, 110)).astype(int))
    hh, cc = np.meshgrid(ha, cells_half, indexing="xy")
    cell_size = args.half_width / cc
    cells_layer = cells_per_hartmann_layer(cell_size, args.half_width, hh)
    labels = classify_hartmann_resolution(
        cell_size,
        args.half_width,
        hh,
        target_cells=args.target_cells,
    )
    alpha = hartmann_channel_effective_drag(args.half_width, hh, args.nu)
    label_code = np.where(labels == "resolved", 0, np.where(labels == "marginal", 1, 2))
    return {
        "ha": ha,
        "cells_half": cells_half,
        "hh": hh,
        "cc": cc,
        "cell_size": cell_size,
        "cells_layer": cells_layer,
        "labels": labels,
        "label_code": label_code,
        "alpha": alpha,
    }


def write_hartmann_csv(outdir: Path, args: argparse.Namespace, scan: dict[str, np.ndarray]) -> None:
    with (outdir / "hartmann_drag_closure_scan.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "Ha",
                "cells_per_half_width",
                "cell_size_over_a",
                "delta_H_over_a",
                "cells_per_Hartmann_layer",
                "resolution_label",
                "alpha_eff_a2_over_nu",
                "closure_recommended",
            ]
        )
        for i, cells in enumerate(scan["cells_half"]):
            for j, ha in enumerate(scan["ha"]):
                label = str(scan["labels"][i, j])
                writer.writerow(
                    [
                        f"{ha:.8e}",
                        f"{cells:d}",
                        f"{1.0 / cells:.8e}",
                        f"{1.0 / ha:.8e}",
                        f"{scan['cells_layer'][i, j]:.8e}",
                        label,
                        f"{scan['alpha'][i, j] * args.half_width**2 / args.nu:.8e}",
                        str(label != "resolved"),
                    ]
                )


def plot_resolution_map(outdir: Path, args: argparse.Namespace, scan: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 5.2))
    mesh = ax.pcolormesh(
        scan["ha"],
        scan["cells_half"],
        scan["cells_layer"],
        norm=LogNorm(vmin=0.1, vmax=max(100.0, float(scan["cells_layer"].max()))),
        shading="auto",
        cmap="viridis",
    )
    ax.contour(scan["ha"], scan["cells_half"], scan["cells_layer"], levels=[4.0, args.target_cells], colors=["w", "k"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Hartmann number Ha")
    ax.set_ylabel("uniform cells per half-width")
    ax.set_title("Hartmann-layer resolution: cells across delta_H")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("cells per Hartmann layer")
    fig.tight_layout()
    fig.savefig(outdir / "hartmann_layer_resolution_map.png", dpi=180)
    plt.close(fig)


def plot_effective_drag(outdir: Path, args: argparse.Namespace, scan: dict[str, np.ndarray]) -> None:
    ha = scan["ha"]
    alpha_norm = hartmann_channel_effective_drag(args.half_width, ha, args.nu) * args.half_width**2 / args.nu
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.loglog(ha, alpha_norm, color="k", lw=1.8, label=r"$\alpha_{\rm eff} a^2/\nu$")
    ax.loglog(ha, ha**2, "--", color="0.45", label=r"$\mathrm{Ha}^2$ large-Ha scale")
    ax.axhline(3.0, color="C3", ls=":", lw=1.1, label="Poiseuille limit 3")
    ax.set_xlabel("Hartmann number Ha")
    ax.set_ylabel(r"normalized effective drag $\alpha_{\rm eff} a^2/\nu$")
    ax.set_title("Mean-flow equivalent drag for insulating Hartmann channel")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, which="both", alpha=0.28)
    fig.tight_layout()
    fig.savefig(outdir / "hartmann_effective_drag_alpha.png", dpi=180)
    plt.close(fig)


def plot_validity_matrix(outdir: Path, scan: dict[str, np.ndarray]) -> None:
    cmap = ListedColormap(["#2b8cbe", "#fdae61", "#d7191c"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(figsize=(7.1, 5.2))
    mesh = ax.pcolormesh(
        scan["ha"],
        scan["cells_half"],
        scan["label_code"],
        cmap=cmap,
        norm=norm,
        shading="auto",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Hartmann number Ha")
    ax.set_ylabel("uniform cells per half-width")
    ax.set_title("Reduced-closure need from layer resolution")
    cbar = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["resolved", "marginal", "closure needed"])
    fig.tight_layout()
    fig.savefig(outdir / "drag_closure_validity_matrix.png", dpi=180)
    plt.close(fig)


def manufactured_feedback(outdir: Path, args: argparse.Namespace) -> dict[str, float | str]:
    r = np.linspace(1.1, 1.9, 85)
    z = np.linspace(-0.45, 0.45, 91)
    rr, zz = np.meshgrid(r, z, indexing="xy")
    j0 = 2.5e6
    jphi = j0 * np.exp(-((rr - 1.5) / 0.18) ** 2 - (zz / 0.16) ** 2)
    jphi -= 0.35 * j0 * np.exp(-((rr - 1.34) / 0.11) ** 2 - ((zz + 0.2) / 0.12) ** 2)

    points = np.column_stack([rr.ravel(), zz.ravel()])
    currents = np.zeros((points.shape[0], 3), dtype=np.float64)
    currents[:, 1] = jphi.ravel()
    snap = CouplingSnapshot(
        points_rz=points,
        current_density=currents,
        time=0.0,
        metadata={"source": "manufactured_current_feedback_diagnostic"},
    )
    h5_path = outdir / "manufactured_current_feedback.h5"
    write_coupling_snapshot(str(h5_path), snap)
    reread = read_coupling_snapshot(str(h5_path))
    jphi_read = reread.toroidal_current_density().reshape(rr.shape)
    feedback = feedback_ratio(sheet_delta_b(jphi_read, args.thickness), args.b_ref)
    label = classify_coupling(rm=0.01, feedback=float(np.max(feedback)), interaction=10.0)
    label_map = np.where(feedback >= 1.0e-3, 1, 0)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0), sharex=True, sharey=True)
    im0 = axes[0].pcolormesh(rr, zz, jphi_read, shading="auto", cmap="coolwarm")
    axes[0].set_title(r"manufactured $J_\phi$ [A/m$^2$]")
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].pcolormesh(rr, zz, feedback, shading="auto", cmap="magma", norm=LogNorm(vmin=1.0e-5, vmax=max(1.0e-3, float(feedback.max()))))
    axes[1].set_title(r"sheet estimate $|\delta B|/|B_{\rm ref}|$")
    fig.colorbar(im1, ax=axes[1])
    im2 = axes[2].pcolormesh(rr, zz, label_map, shading="auto", cmap=ListedColormap(["#4575b4", "#d73027"]), vmin=0, vmax=1)
    axes[2].set_title("feedback flag >= 1e-3")
    cbar = fig.colorbar(im2, ax=axes[2], ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(["below", "candidate"])
    for ax in axes:
        ax.set_xlabel("R [m]")
        ax.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("Z [m]")
    fig.suptitle(f"CouplingSnapshot HDF5 roundtrip; max feedback={feedback.max():.2e}; label={label}")
    fig.tight_layout()
    fig.savefig(outdir / "manufactured_current_feedback.png", dpi=180)
    plt.close(fig)

    return {
        "h5": str(h5_path),
        "max_abs_jphi": float(np.max(np.abs(jphi_read))),
        "max_feedback": float(np.max(feedback)),
        "classification": label,
    }


def write_summary(outdir: Path, args: argparse.Namespace, scan: dict[str, np.ndarray], feedback_info: dict[str, float | str]) -> None:
    labels, counts = np.unique(scan["labels"], return_counts=True)
    label_counts = {str(label): int(count) for label, count in zip(labels, counts)}
    lines = [
        "# Hartmann effective-drag closure diagnostics",
        "",
        "This is an analytic/reduced-model infrastructure artifact. It is not a MUG run,",
        "not a TokaMaker coupling result, and not a claim of standalone solver validation.",
        "",
        "Purpose:",
        "- Identify when a full-device mesh under-resolves delta_H = a/Ha.",
        "- Provide the mean-flow equivalent alpha_eff for an insulating Hartmann channel.",
        "- Exercise the CouplingSnapshot HDF5 path with a manufactured current distribution.",
        "",
        "Assumptions:",
        f"- Channel half-width a = {args.half_width:g} m.",
        f"- Kinematic viscosity nu = {args.nu:g} m^2/s.",
        f"- Resolved threshold = {args.target_cells:g} cells per Hartmann layer.",
        "- Marginal means half the target threshold to the target threshold.",
        "- Reduced drag is intended for full-device underresolved-layer meshes, not for",
        "  replacing resolved standalone Hartmann verification.",
        "",
        "Resolution matrix counts:",
        "",
        "| label | count |",
        "|---|---:|",
    ]
    for label in ["resolved", "marginal", "underresolved"]:
        lines.append(f"| {label} | {label_counts.get(label, 0)} |")
    lines += [
        "",
        "Manufactured current-feedback diagnostic:",
        "",
        f"- HDF5 file: `{Path(str(feedback_info['h5'])).name}`",
        f"- max |J_phi| = {float(feedback_info['max_abs_jphi']):.3e} A/m^2",
        f"- max sheet-estimate |delta B|/|B_ref| = {float(feedback_info['max_feedback']):.3e}",
        f"- triage classification = `{feedback_info['classification']}`",
        "",
        "Generated plots:",
        "- `hartmann_layer_resolution_map.png`",
        "- `hartmann_effective_drag_alpha.png`",
        "- `drag_closure_validity_matrix.png`",
        "- `manufactured_current_feedback.png`",
        "",
        "Red-team notes:",
        "- alpha_eff is calibrated to the mean velocity/forcing balance only.",
        "- It must not be used to claim pointwise boundary-layer profile agreement.",
        "- A feedback flag from the manufactured current is only a schema/triage check.",
    ]
    (outdir / "hartmann_drag_closure_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    scan = build_hartmann_scan(args)
    write_hartmann_csv(outdir, args, scan)
    plot_resolution_map(outdir, args, scan)
    plot_effective_drag(outdir, args, scan)
    plot_validity_matrix(outdir, scan)
    feedback_info = manufactured_feedback(outdir, args)
    write_summary(outdir, args, scan, feedback_info)
    print(f"Wrote Hartmann closure diagnostics to {outdir}")


if __name__ == "__main__":
    main()
