#!/usr/bin/env python3
"""Homogeneous-Neumann electric-potential Poisson MMS diagnostic.

This verifies a mean-gauged finite-difference solve for laplacian(phi)=source
with zero normal derivative on all rectangle walls. It is closer to insulating
electric-potential wall checks than the Dirichlet MMS, but still not a MUG or
finite-element result.
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
    return parser.parse_args()


def exact_solution(x_grid: np.ndarray, z_grid: np.ndarray) -> np.ndarray:
    return np.cos(np.pi * x_grid) * np.cos(np.pi * z_grid)


def exact_source(x_grid: np.ndarray, z_grid: np.ndarray) -> np.ndarray:
    return -2.0 * np.pi**2 * exact_solution(x_grid, z_grid)


def solve_case(n: int) -> dict[str, np.ndarray | float]:
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.0, 1.0, n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = exact_solution(x_grid, z_grid)
    phi = solve_potential_neumann_2d(exact_source(x_grid, z_grid), x, z, mean_value=0.0)
    error = phi - phi_exact
    error -= np.mean(error)
    return {
        "n": float(n),
        "h": float(x[1] - x[0]),
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "phi_exact": phi_exact,
        "phi": phi,
        "error": error,
        "max_error": float(np.max(np.abs(error))),
        "rms_error": float(np.sqrt(np.mean(error**2))),
        "mean_phi": float(np.mean(phi)),
    }


def run_cases() -> list[dict[str, np.ndarray | float]]:
    return [solve_case(n) for n in (17, 25, 33, 49, 65)]


def write_csv(outdir: Path, cases: list[dict[str, np.ndarray | float]]) -> None:
    path = outdir / "inductionless_neumann_mms_metrics.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["n", "h", "max_error", "rms_error", "mean_phi", "observed_order_max"])
        previous = None
        for case in cases:
            if previous is None:
                order_max = ""
            else:
                order_max_value = np.log(float(previous["max_error"]) / float(case["max_error"])) / np.log(
                    float(previous["h"]) / float(case["h"])
                )
                order_max = f"{order_max_value:.8e}"
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


def plot_solution(outdir: Path, case: dict[str, np.ndarray | float]) -> None:
    x_grid = case["x_grid"]
    z_grid = case["z_grid"]
    phi_exact = case["phi_exact"]
    phi = case["phi"]
    error = case["error"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharex=True, sharey=True)
    panels = [
        (phi_exact, r"exact $\phi$", "viridis"),
        (phi, r"solved $\phi$", "viridis"),
        (error, r"mean-free error", "coolwarm"),
    ]
    for ax, (values, title, cmap) in zip(axes, panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    fig.suptitle(f"Homogeneous Neumann Poisson MMS, n={int(case['n'])}")
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_neumann_mms_solution.png", dpi=180)
    plt.close(fig)


def plot_convergence(outdir: Path, cases: list[dict[str, np.ndarray | float]]) -> None:
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
    ax.set_ylabel(r"mean-free $\phi$ error")
    ax.set_title("Homogeneous Neumann Poisson MMS convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "inductionless_neumann_mms_convergence.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, cases: list[dict[str, np.ndarray | float]]) -> None:
    finest = cases[-1]
    lines = [
        "# Homogeneous Neumann electric-potential Poisson MMS",
        "",
        "This verifies a mean-gauged scalar finite-difference reference solve with",
        "zero normal derivative on all rectangle walls. It is not a MUG result.",
        "",
        "Equation and manufactured solution:",
        "- Solve `laplacian(phi) = source` on `[0,1] x [0,1]`.",
        "- `phi = cos(pi x) cos(pi z)`, which has homogeneous Neumann walls.",
        "- `source = -2 pi^2 phi`, whose grid mean is zero to roundoff.",
        "- The nullspace is fixed by enforcing mean(phi) = 0.",
        "",
        "Finest-grid result:",
        f"- n = {int(finest['n'])}.",
        f"- h = {float(finest['h']):.3e}.",
        f"- max mean-free error = {float(finest['max_error']):.3e}.",
        f"- rms mean-free error = {float(finest['rms_error']):.3e}.",
        f"- solved mean(phi) = {float(finest['mean_phi']):.3e}.",
        "",
        "Generated files:",
        "- `inductionless_neumann_mms_metrics.csv`",
        "- `inductionless_neumann_mms_solution.png`",
        "- `inductionless_neumann_mms_convergence.png`",
        "",
        "Red-team interpretation:",
        "- This is homogeneous Neumann only: no nonzero motional-wall Neumann data.",
        "- It has no Robin conducting-wall closure and no material interfaces.",
        "- It is a finite-difference reference, not the finite-element equation needed in MUG.",
    ]
    (outdir / "inductionless_neumann_mms_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cases = run_cases()
    write_csv(outdir, cases)
    plot_solution(outdir, cases[-1])
    plot_convergence(outdir, cases)
    write_summary(outdir, cases)
    print(f"Wrote inductionless Neumann MMS diagnostics to {outdir}")


if __name__ == "__main__":
    main()
