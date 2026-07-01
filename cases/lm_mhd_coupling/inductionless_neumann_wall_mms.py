#!/usr/bin/env python3
"""Nonzero-Neumann electric-potential wall-data MMS diagnostic.

This verifies the reference solve for laplacian(phi)=source with prescribed
outward-normal derivative on each rectangle wall. It is a scalar MMS for the
kind of boundary data used by insulating inductionless walls,
dphi/dn = (u x B).n, not a full LM-MHD solve.
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
    from OpenFUSIONToolkit.LM_MHD import solve_potential_neumann_2d  # noqa: E402
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    solve_potential_neumann_2d = inductionless.solve_potential_neumann_2d


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--n", type=int, default=65, help="grid points in each direction")
    return parser.parse_args()


def build_case(n: int) -> dict[str, np.ndarray | float | dict[str, np.ndarray]]:
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.0, 1.0, n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = (x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2
    phi_exact -= np.mean(phi_exact)
    source = np.full_like(phi_exact, 4.0)
    normal_gradient = {
        "x_min": np.ones_like(z),
        "x_max": np.ones_like(z),
        "z_min": np.ones_like(x),
        "z_max": np.ones_like(x),
    }
    phi = solve_potential_neumann_2d(source, x, z, mean_value=0.0, normal_gradient=normal_gradient)
    error = phi - phi_exact
    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "phi_exact": phi_exact,
        "phi": phi,
        "error": error,
        "source": source,
        "normal_gradient": normal_gradient,
        "max_error": float(np.max(np.abs(error))),
        "rms_error": float(np.sqrt(np.mean(error**2))),
        "mean_phi": float(np.mean(phi)),
    }


def write_metrics(outdir: Path, case: dict[str, np.ndarray | float | dict[str, np.ndarray]]) -> None:
    path = outdir / "inductionless_neumann_wall_mms_metrics.csv"
    rows = [
        ("n", int(len(case["x"]))),
        ("h", float(case["x"][1] - case["x"][0])),
        ("source", 4.0),
        ("wall_normal_gradient", 1.0),
        ("max_error", float(case["max_error"])),
        ("rms_error", float(case["rms_error"])),
        ("mean_phi", float(case["mean_phi"])),
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in rows:
            writer.writerow([key, f"{value:.12e}" if isinstance(value, float) else value])


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float | dict[str, np.ndarray]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharex=True, sharey=True)
    panels = [
        (case["phi_exact"], r"exact $\phi$", "viridis"),
        (case["phi"], r"solved $\phi$", "viridis"),
        (case["error"], r"error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Nonzero-Neumann wall-data Poisson MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_neumann_wall_mms_solution.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, case: dict[str, np.ndarray | float | dict[str, np.ndarray]]) -> None:
    lines = [
        "# Nonzero-Neumann wall-data Poisson MMS",
        "",
        "This verifies prescribed outward-normal derivative data for the scalar",
        "electric-potential reference solve. It is not a MUG result.",
        "",
        "Manufactured solution:",
        "- `phi = (x - 1/2)^2 + (z - 1/2)^2 - mean(phi)`.",
        "- `laplacian(phi) = 4`.",
        "- outward `dphi/dn = 1` on all four rectangle walls.",
        "- mean(phi) is fixed to zero as the Neumann gauge.",
        "",
        "Result:",
        f"- n = {int(len(case['x']))}.",
        f"- max error = {float(case['max_error']):.3e}.",
        f"- rms error = {float(case['rms_error']):.3e}.",
        f"- solved mean(phi) = {float(case['mean_phi']):.3e}.",
        "",
        "Red-team interpretation:",
        "- This checks scalar wall-gradient plumbing only.",
        "- It does not yet set the wall gradient from an actual `u x B` field.",
        "- It does not include Robin wall conductance or material interfaces.",
    ]
    (outdir / "inductionless_neumann_wall_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    case = build_case(args.n)
    write_metrics(outdir, case)
    plot_solution(outdir, case)
    write_summary(outdir, case)
    print(f"Wrote nonzero-Neumann wall MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
