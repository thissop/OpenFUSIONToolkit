#!/usr/bin/env python3
"""Two-region conductivity-jump potential MMS diagnostic.

This verifies that the Python variable-conductivity reference can represent a
face-aligned material jump with continuity of potential and normal current. It
is a finite-difference reference artifact, not a MUG result.
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
    from OpenFUSIONToolkit.LM_MHD import solve_variable_conductivity_dirichlet_2d  # noqa: E402
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    solve_variable_conductivity_dirichlet_2d = inductionless.solve_variable_conductivity_dirichlet_2d


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--nx", type=int, default=64, help="even number of x grid points")
    parser.add_argument("--nz", type=int, default=33, help="number of z grid points")
    parser.add_argument("--sigma-left", type=float, default=1.0, help="left-region conductivity")
    parser.add_argument("--sigma-right", type=float, default=4.0, help="right-region conductivity")
    return parser.parse_args()


def build_case(args: argparse.Namespace) -> dict[str, np.ndarray | float]:
    if args.nx % 2 != 0:
        raise ValueError("--nx must be even so the interface lies on a grid face")
    x = np.linspace(0.0, 1.0, args.nx)
    z = np.linspace(0.0, 1.0, args.nz)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    interface = 0.5
    sigma_left = args.sigma_left
    sigma_right = args.sigma_right
    flux = 1.0 / (interface / sigma_left + (1.0 - interface) / sigma_right)
    conductivity = np.where(x_grid < interface, sigma_left, sigma_right)
    phi_1d = np.where(
        x <= interface,
        flux * x / sigma_left,
        flux * interface / sigma_left + flux * (x - interface) / sigma_right,
    )
    phi_exact = np.broadcast_to(phi_1d, x_grid.shape)
    source = np.zeros_like(phi_exact)
    phi = solve_variable_conductivity_dirichlet_2d(conductivity, source, x, z, boundary_values=phi_exact)
    error = phi - phi_exact
    left_idx = x.size // 2 - 1
    right_idx = x.size // 2
    h = x[1] - x[0]
    sigma_face = 2.0 * sigma_left * sigma_right / (sigma_left + sigma_right)
    interface_flux = sigma_face * (phi[:, right_idx] - phi[:, left_idx]) / h
    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "conductivity": conductivity,
        "phi_exact": phi_exact,
        "phi": phi,
        "error": error,
        "interface_flux": interface_flux,
        "target_flux": float(flux),
        "max_error": float(np.max(np.abs(error))),
        "rms_error": float(np.sqrt(np.mean(error**2))),
        "max_flux_error": float(np.max(np.abs(interface_flux - flux))),
        "interface": interface,
    }


def write_metrics(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    rows = [
        ("target_flux", float(case["target_flux"])),
        ("max_error", float(case["max_error"])),
        ("rms_error", float(case["rms_error"])),
        ("max_interface_flux_error", float(case["max_flux_error"])),
    ]
    with (outdir / "inductionless_two_region_sigma_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in rows:
            writer.writerow([key, f"{value:.12e}"])


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    panels = [
        (case["conductivity"], r"$\sigma$", "viridis"),
        (case["phi_exact"], r"exact $\phi$", "magma"),
        (case["phi"], r"solved $\phi$", "magma"),
        (case["error"], r"$\phi$ error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.axvline(float(case["interface"]), color="white", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle("Two-region conductivity-jump MMS")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_two_region_sigma_mms_solution.png", dpi=180)
    plt.close(fig)


def plot_flux(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    fig, ax = plt.subplots(figsize=(6.3, 4.2))
    z = case["z"]
    flux = case["interface_flux"]
    target = float(case["target_flux"])
    ax.plot(z, flux, label="interface flux")
    ax.axhline(target, color="0.35", ls="--", label="target flux")
    ax.set_xlabel("z")
    ax.set_ylabel(r"$\sigma \partial\phi/\partial x$")
    ax.set_title("Normal current continuity at material interface")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_two_region_sigma_mms_flux.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    lines = [
        "# Two-region conductivity-jump potential MMS",
        "",
        "This is a finite-difference reference diagnostic for material-interface",
        "coefficient handling. It does not run MUG.",
        "",
        "Manufactured setup:",
        "- Two conductivities: sigma_left = 1, sigma_right = 4.",
        "- The interface is aligned with an x-grid face at x = 0.5.",
        "- The exact potential is piecewise linear.",
        "- The source is zero.",
        "- Potential and normal current are continuous across the interface.",
        "",
        "Result:",
        f"- target interface flux = {float(case['target_flux']):.6e}.",
        f"- max potential error = {float(case['max_error']):.3e}.",
        f"- max interface-flux error = {float(case['max_flux_error']):.3e}.",
        "",
        "Generated files:",
        "- `inductionless_two_region_sigma_mms_metrics.csv`",
        "- `inductionless_two_region_sigma_mms_solution.png`",
        "- `inductionless_two_region_sigma_mms_flux.png`",
        "",
        "Red-team interpretation:",
        "- The interface is face-aligned and one-dimensional in x.",
        "- This checks harmonic face averaging and flux continuity only.",
        "- It is not a finite-element multi-region implementation.",
    ]
    (outdir / "inductionless_two_region_sigma_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    case = build_case(args)
    write_metrics(outdir, case)
    plot_solution(outdir, case)
    plot_flux(outdir, case)
    write_summary(outdir, case)
    print(f"Wrote two-region conductivity MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
