#!/usr/bin/env python3
"""Offline inductionless current diagnostics for structured CSV fields.

This is a bridge utility for future MUG/TokaMaker postprocessing: read a
structured x-z CSV containing velocity, magnetic field, conductivity, and
potential; reconstruct `J = sigma(-grad(phi) + u x B)`; then report `div J`,
wall-normal current, and optional target-current errors. It is not a solver.
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
        charge_conservation_residual_2d,
        divergence_free_current_from_streamfunction,
        reconstruct_inductionless_current_2d,
        structured_gradient_2d,
        wall_normal_current_extrema,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    inductionless_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    charge_conservation_residual_2d = inductionless.charge_conservation_residual_2d
    divergence_free_current_from_streamfunction = inductionless.divergence_free_current_from_streamfunction
    reconstruct_inductionless_current_2d = inductionless.reconstruct_inductionless_current_2d
    structured_gradient_2d = inductionless.structured_gradient_2d
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", nargs="?", help="input CSV; omit with --demo")
    parser.add_argument("--demo", action="store_true", help="generate and postprocess a manufactured CSV")
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--prefix", default="inductionless_postprocess_demo", help="output filename prefix")
    parser.add_argument("--n", type=int, default=65, help="demo grid points in each direction")
    parser.add_argument("--b-y", type=float, default=4.0, help="fallback/demo transverse magnetic field")
    return parser.parse_args()


def demo_rows(n: int, b_y: float) -> dict[str, np.ndarray]:
    x = np.linspace(0.0, 1.0, n)
    z = np.linspace(0.0, 1.0, n)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    conductivity = 1.0 + 0.25 * x_grid + 0.15 * z_grid
    phi = 0.015 * ((x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2)
    phi -= np.mean(phi)
    streamfunction = 0.02 * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)
    emf = structured_gradient_2d(phi, x, z) + target_current / conductivity[..., None]
    velocity = np.zeros_like(target_current)
    velocity[..., 0] = emf[..., 2] / b_y
    velocity[..., 2] = -emf[..., 0] / b_y
    return {
        "x": x_grid.ravel(),
        "z": z_grid.ravel(),
        "u_x": velocity[..., 0].ravel(),
        "u_y": velocity[..., 1].ravel(),
        "u_z": velocity[..., 2].ravel(),
        "b_x": np.zeros(x_grid.size),
        "b_y": np.full(x_grid.size, b_y),
        "b_z": np.zeros(x_grid.size),
        "sigma": conductivity.ravel(),
        "phi": phi.ravel(),
        "target_jx": target_current[..., 0].ravel(),
        "target_jy": target_current[..., 1].ravel(),
        "target_jz": target_current[..., 2].ravel(),
    }


def write_demo_csv(path: Path, rows: dict[str, np.ndarray]) -> None:
    names = list(rows)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(names)
        for values in zip(*(rows[name] for name in names)):
            writer.writerow([f"{float(value):.16e}" for value in values])


def load_structured_csv(path: Path, fallback_b_y: float) -> dict[str, np.ndarray | None]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.dtype.names is None:
        raise ValueError("CSV must have a header row")
    names = set(data.dtype.names)
    required = {"x", "z", "u_x", "u_z", "sigma", "phi"}
    missing = required - names
    if missing:
        raise ValueError(f"CSV missing required column(s): {sorted(missing)}")
    x = np.unique(data["x"])
    z = np.unique(data["z"])
    if x.size * z.size != data.shape[0]:
        raise ValueError("CSV rows must form a complete tensor-product x-z grid")
    order = np.lexsort((data["x"], data["z"]))
    sorted_data = data[order]
    x_grid = sorted_data["x"].reshape(z.size, x.size)
    z_grid = sorted_data["z"].reshape(z.size, x.size)
    if not np.allclose(x_grid, x[None, :]) or not np.allclose(z_grid, z[:, None]):
        raise ValueError("CSV rows do not form a structured x-z grid after sorting")

    velocity = np.zeros((z.size, x.size, 3), dtype=np.float64)
    velocity[..., 0] = sorted_data["u_x"].reshape(z.size, x.size)
    velocity[..., 1] = _optional_grid(sorted_data, "u_y", z.size, x.size)
    velocity[..., 2] = sorted_data["u_z"].reshape(z.size, x.size)

    magnetic_field = np.zeros_like(velocity)
    magnetic_field[..., 0] = _optional_grid(sorted_data, "b_x", z.size, x.size)
    magnetic_field[..., 1] = _optional_grid(sorted_data, "b_y", z.size, x.size, fallback=fallback_b_y)
    magnetic_field[..., 2] = _optional_grid(sorted_data, "b_z", z.size, x.size)

    target_current = None
    if {"target_jx", "target_jz"} <= names:
        target_current = np.zeros_like(velocity)
        target_current[..., 0] = sorted_data["target_jx"].reshape(z.size, x.size)
        target_current[..., 1] = _optional_grid(sorted_data, "target_jy", z.size, x.size)
        target_current[..., 2] = sorted_data["target_jz"].reshape(z.size, x.size)

    return {
        "x": x,
        "z": z,
        "x_grid": x_grid,
        "z_grid": z_grid,
        "velocity": velocity,
        "magnetic_field": magnetic_field,
        "conductivity": sorted_data["sigma"].reshape(z.size, x.size),
        "potential": sorted_data["phi"].reshape(z.size, x.size),
        "target_current": target_current,
    }


def _optional_grid(
    data: np.ndarray,
    name: str,
    nz: int,
    nx: int,
    fallback: float = 0.0,
) -> np.ndarray:
    if name in data.dtype.names:
        return data[name].reshape(nz, nx)
    return np.full((nz, nx), fallback, dtype=np.float64)


def postprocess(fields: dict[str, np.ndarray | None]) -> dict[str, np.ndarray | dict[str, float] | None]:
    current = reconstruct_inductionless_current_2d(
        fields["potential"],
        fields["velocity"],
        fields["magnetic_field"],
        fields["conductivity"],
        fields["x"],
        fields["z"],
    )
    residual = charge_conservation_residual_2d(current, fields["x"], fields["z"])
    wall_current = wall_normal_current_extrema(current)
    target = fields["target_current"]
    error = None if target is None else current - target
    return {
        **fields,
        "current": current,
        "residual": residual,
        "wall_current": wall_current,
        "current_error": error,
    }


def metric_values(result: dict[str, np.ndarray | dict[str, float] | None]) -> dict[str, float]:
    current = result["current"]
    residual = result["residual"]
    wall_current = result["wall_current"]
    assert isinstance(wall_current, dict)
    interior = residual[2:-2, 2:-2] if residual.shape[0] > 4 and residual.shape[1] > 4 else residual
    metrics = {
        "max_abs_current": float(np.max(np.linalg.norm(current, axis=-1))),
        "max_abs_divJ_interior": float(np.max(np.abs(interior))),
        "rms_divJ_interior": float(np.sqrt(np.mean(interior**2))),
        "max_wall_normal_current": float(max(wall_current.values())),
    }
    if result["current_error"] is not None:
        error = result["current_error"]
        metrics["max_abs_current_error"] = float(np.max(np.abs(error)))
        metrics["rms_current_error"] = float(np.sqrt(np.mean(error**2)))
    return metrics


def write_metrics(outdir: Path, prefix: str, metrics: dict[str, float]) -> None:
    with (outdir / f"{prefix}_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, f"{value:.12e}"])


def plot_result(outdir: Path, prefix: str, result: dict[str, np.ndarray | dict[str, float] | None]) -> None:
    x_grid = result["x_grid"]
    z_grid = result["z_grid"]
    current = result["current"]
    residual = result["residual"]
    panels = [
        (result["conductivity"], r"$\sigma$", "viridis"),
        (result["potential"], r"$\phi$", "magma"),
        (np.linalg.norm(current, axis=-1), r"$|J|$", "viridis"),
        (residual, r"$\nabla\cdot J$", "coolwarm"),
    ]
    if result["current_error"] is not None:
        panels.extend(
            [
                (np.linalg.norm(result["target_current"], axis=-1), r"target $|J|$", "viridis"),
                (np.linalg.norm(result["current_error"], axis=-1), r"$|J-J_{target}|$", "magma"),
            ]
        )
    ncols = 3 if len(panels) > 4 else 2
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 3.7 * nrows), sharex=True, sharey=True)
    for ax, (values, title, cmap) in zip(np.ravel(axes), panels):
        img = ax.pcolormesh(x_grid, z_grid, values, shading="auto", cmap=cmap)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_title(title)
        fig.colorbar(img, ax=ax)
    for ax in np.ravel(axes)[len(panels) :]:
        ax.axis("off")
    fig.suptitle("Offline inductionless CSV diagnostics")
    fig.tight_layout()
    fig.savefig(outdir / f"{prefix}_current.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, prefix: str, metrics: dict[str, float], input_path: Path) -> None:
    lines = [
        "# Offline inductionless CSV diagnostics",
        "",
        f"Input CSV: `{input_path}`",
        "",
        "This postprocesses structured x-z fields into inductionless current",
        "diagnostics. It does not solve for the potential and does not run MUG.",
        "",
        "Metrics:",
    ]
    lines.extend(f"- {key} = {value:.3e}." for key, value in metrics.items())
    lines.extend(
        [
            "",
            "Expected CSV columns:",
            "- Required: `x`, `z`, `u_x`, `u_z`, `sigma`, `phi`.",
            "- Optional: `u_y`, `b_x`, `b_y`, `b_z`, `target_jx`, `target_jy`, `target_jz`.",
            "",
            "Red-team interpretation:",
            "- Structured uniform-grid diagnostic only.",
            "- Not a finite-element projection and not a conservative interface-current",
            "  operator at material jumps.",
            "- For MUG use, first export/interpolate fields onto a documented structured",
            "  grid and keep the interpolation script with the result.",
        ]
    )
    (outdir / f"{prefix}_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    if args.demo:
        rows = demo_rows(args.n, args.b_y)
        input_path = outdir / f"{args.prefix}_input.csv"
        write_demo_csv(input_path, rows)
    elif args.csv_file:
        input_path = Path(args.csv_file).resolve()
    else:
        raise SystemExit("Provide a CSV file or use --demo")

    fields = load_structured_csv(input_path, fallback_b_y=args.b_y)
    result = postprocess(fields)
    metrics = metric_values(result)
    write_metrics(outdir, args.prefix, metrics)
    plot_result(outdir, args.prefix, result)
    write_summary(outdir, args.prefix, metrics, input_path)
    print(f"Wrote offline inductionless CSV diagnostics to {outdir}")


if __name__ == "__main__":
    main()
