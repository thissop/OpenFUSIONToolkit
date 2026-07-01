#!/usr/bin/env python3
"""Variable-conductivity Neumann-flux potential MMS diagnostic.

This verifies a smooth variable-conductivity reference solve for
div(sigma grad(phi)) = source with prescribed outward normal flux
sigma*dphi/dn on all walls. It is a Python reference artifact, not a MUG solve.
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
    from OpenFUSIONToolkit.LM_MHD import solve_variable_conductivity_neumann_2d  # noqa: E402
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    solve_variable_conductivity_neumann_2d = inductionless.solve_variable_conductivity_neumann_2d


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    return parser.parse_args()


def build_case(n: int) -> dict[str, np.ndarray | float | dict[str, np.ndarray]]:
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.0, 1.0, n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = (x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2
    phi_exact -= np.mean(phi_exact)
    conductivity = 1.0 + 0.3 * x_grid + 0.2 * z_grid
    source = 4.0 * conductivity + 0.6 * (x_grid - 0.5) + 0.4 * (z_grid - 0.5)
    normal_flux = {
        "x_min": conductivity[:, 0],
        "x_max": conductivity[:, -1],
        "z_min": conductivity[0, :],
        "z_max": conductivity[-1, :],
    }
    phi = solve_variable_conductivity_neumann_2d(
        conductivity,
        source,
        x,
        z,
        mean_value=0.0,
        normal_flux=normal_flux,
    )
    error = phi - phi_exact
    return {
        "n": float(n),
        "h": float(x[1] - x[0]),
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "conductivity": conductivity,
        "source": source,
        "normal_flux": normal_flux,
        "phi_exact": phi_exact,
        "phi": phi,
        "error": error,
        "max_error": float(np.max(np.abs(error))),
        "rms_error": float(np.sqrt(np.mean(error**2))),
        "mean_phi": float(np.mean(phi)),
    }


def run_cases() -> list[dict[str, np.ndarray | float | dict[str, np.ndarray]]]:
    return [build_case(n) for n in (17, 25, 33, 49, 65)]


def write_csv(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, np.ndarray]]]) -> None:
    with (outdir / "inductionless_variable_sigma_neumann_mms_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["n", "h", "max_error", "rms_error", "mean_phi", "observed_order_max"])
        previous = None
        for case in cases:
            if previous is None:
                order_max = ""
            else:
                order = np.log(float(previous["max_error"]) / float(case["max_error"])) / np.log(
                    float(previous["h"]) / float(case["h"])
                )
                order_max = f"{order:.8e}"
            writer.writerow(
                [
                    int(case["n"]),
                    f"{float(case['h']):.8e}",
                    f"{float(case['max_error']):.8e}",
                    f"{float(case['rms_error']):.8e}",
                    f"{float(case['mean_phi']):.8e}",
                    order_max,
                ]
            )
            previous = case


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float | dict[str, np.ndarray]]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    fig, axes = plt.subplots(1, 4, figsize=(14.8, 3.8), sharex=True, sharey=True)
    panels = [
        (case["conductivity"], r"$\sigma(x,z)$", "viridis"),
        (case["phi_exact"], r"exact $\phi$", "magma"),
        (case["phi"], r"solved $\phi$", "magma"),
        (case["error"], r"$\phi$ error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle(f"Variable-conductivity Neumann-flux MMS, n={int(case['n'])}")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_variable_sigma_neumann_mms_solution.png", dpi=180)
    plt.close(fig)


def plot_convergence(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, np.ndarray]]]) -> None:
    h = np.array([float(case["h"]) for case in cases])
    max_error = np.array([float(case["max_error"]) for case in cases])
    rms_error = np.array([float(case["rms_error"]) for case in cases])
    reference = max_error[0] * (h / h[0]) ** 2

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(h, max_error, "o-", label="max error")
    ax.loglog(h, rms_error, "s-", label="rms error")
    ax.loglog(h, reference, "--", color="0.4", label=r"$O(h^2)$ reference")
    ax.invert_xaxis()
    ax.set_xlabel("grid spacing h")
    ax.set_ylabel(r"mean-gauged $\phi$ error")
    ax.set_title("Variable-conductivity Neumann MMS convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_variable_sigma_neumann_mms_convergence.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, cases: list[dict[str, np.ndarray | float | dict[str, np.ndarray]]]) -> None:
    finest = cases[-1]
    lines = [
        "# Variable-conductivity Neumann-flux potential MMS",
        "",
        "This is a Python reference diagnostic for smooth variable-conductivity",
        "Neumann flux closure. It does not run MUG.",
        "",
        "Manufactured setup:",
        "- `sigma = 1 + 0.3 x + 0.2 z`.",
        "- `phi = (x - 1/2)^2 + (z - 1/2)^2 - mean(phi)`.",
        "- `source = div(sigma grad(phi))` evaluated analytically.",
        "- outward wall flux is `sigma dphi/dn`.",
        "- mean(phi) is fixed to zero.",
        "",
        "Finest-grid result:",
        f"- n = {int(finest['n'])}.",
        f"- max error = {float(finest['max_error']):.3e}.",
        f"- rms error = {float(finest['rms_error']):.3e}.",
        f"- solved mean(phi) = {float(finest['mean_phi']):.3e}.",
        "",
        "Generated files:",
        "- `inductionless_variable_sigma_neumann_mms_metrics.csv`",
        "- `inductionless_variable_sigma_neumann_mms_solution.png`",
        "- `inductionless_variable_sigma_neumann_mms_convergence.png`",
        "",
        "Red-team interpretation:",
        "- Smooth scalar conductivity only.",
        "- Prescribed flux data only; not generated from `u x B` yet.",
        "- No discontinuous material interface in this Neumann diagnostic.",
        "- Uniform-grid finite differences only, not finite elements.",
    ]
    (outdir / "inductionless_variable_sigma_neumann_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cases = run_cases()
    write_csv(outdir, cases)
    plot_solution(outdir, cases[-1])
    plot_convergence(outdir, cases)
    write_summary(outdir, cases)
    print(f"Wrote variable-conductivity Neumann MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
