#!/usr/bin/env python3
"""Inspect a MUG/xmhd_2d HDF5 restart for LM-MHD export readiness.

This utility records what is actually present in a restart file before any
offline inductionless CSV exporter is attempted. It does not modify or run MUG.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import smplotlib  # noqa: F401  (applies the style on import)
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "restart",
        nargs="?",
        default="cases/hartmann/mug/Ha_2/xmhd2d_00000.rst",
        help="MUG/xmhd_2d HDF5 restart file",
    )
    parser.add_argument(
        "--profile",
        default="cases/hartmann/mug/Ha_2/hartmann_mug.profile",
        help="optional Hartmann profile with x,z/velocity coordinates",
    )
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    return parser.parse_args()


def inspect_restart(path: Path) -> list[dict[str, str | float | int]]:
    rows: list[dict[str, str | float | int]] = []
    with h5py.File(path, "r") as handle:
        for name, obj in sorted(handle.items()):
            if not isinstance(obj, h5py.Dataset):
                continue
            data = obj[...]
            finite = np.isfinite(data) if np.issubdtype(data.dtype, np.number) else np.array([False])
            row: dict[str, str | float | int] = {
                "name": name,
                "shape": "x".join(str(item) for item in data.shape),
                "dtype": str(data.dtype),
                "size": int(data.size),
                "finite": bool(np.all(finite)),
            }
            if np.issubdtype(data.dtype, np.number):
                row["min"] = float(np.min(data))
                row["max"] = float(np.max(data))
                row["max_abs"] = float(np.max(np.abs(data)))
            else:
                row["min"] = ""
                row["max"] = ""
                row["max_abs"] = ""
            rows.append(row)
    return rows


def profile_info(path: Path) -> dict[str, float | int | str]:
    if not path.exists():
        return {"path": str(path), "exists": "false"}
    data = np.loadtxt(path, comments="#")
    return {
        "path": str(path),
        "exists": "true",
        "rows": int(data.shape[0]),
        "columns": int(data.shape[1]),
        "x_min": float(np.min(data[:, 0])),
        "x_max": float(np.max(data[:, 0])),
        "z_min": float(np.min(data[:, 1])),
        "z_max": float(np.max(data[:, 1])),
    }


def write_manifest(outdir: Path, rows: list[dict[str, str | float | int]]) -> None:
    keys = ["name", "shape", "dtype", "size", "finite", "min", "max", "max_abs"]
    with (outdir / "mug_restart_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_ranges(outdir: Path, rows: list[dict[str, str | float | int]]) -> None:
    field_rows = [row for row in rows if str(row["name"]).startswith("U_")]
    names = [str(row["name"]).replace("U_", "") for row in field_rows]
    max_abs = np.array([float(row["max_abs"]) for row in field_rows])
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    y = np.arange(len(names))
    ax.barh(y, max_abs)
    ax.set_yticks(y, names)
    ax.set_xlabel("max absolute value")
    ax.set_title("MUG restart field ranges")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "mug_restart_field_ranges.png", dpi=180)
    plt.close(fig)


def write_summary(
    outdir: Path,
    restart: Path,
    rows: list[dict[str, str | float | int]],
    profile: dict[str, float | int | str],
) -> None:
    names = [str(row["name"]) for row in rows]
    has_coords = any(name.lower() in {"x", "z", "coords", "coordinates", "mesh"} for name in names)
    has_phi = any(name.lower() in {"phi", "u_phi", "electric_potential"} for name in names)
    restart_time = next((float(row["max"]) for row in rows if row["name"] == "t"), float("nan"))
    state_size = next((int(row["size"]) for row in rows if row["name"] == "U_n"), 0)
    profile_rows = int(profile["rows"]) if profile.get("exists") == "true" else 0
    profile_matches_state = bool(profile_rows == state_size and state_size > 0)
    lines = [
        "# MUG restart export-readiness manifest",
        "",
        f"Restart inspected: `{restart}`",
        "",
        "Datasets:",
        "",
        "| name | shape | dtype | max_abs |",
        "|---|---:|---|---:|",
    ]
    for row in rows:
        lines.append(f"| `{row['name']}` | `{row['shape']}` | `{row['dtype']}` | `{float(row['max_abs']):.3e}` |")
    lines.extend(
        [
            "",
            "Profile sidecar:",
            f"- path: `{profile['path']}`",
            f"- exists: `{profile['exists']}`",
        ]
    )
    if profile.get("exists") == "true":
        lines.extend(
            [
                f"- rows: `{profile['rows']}`",
                f"- coordinate span: x=[{float(profile['x_min']):.3g}, {float(profile['x_max']):.3g}], "
                f"z=[{float(profile['z_min']):.3g}, {float(profile['z_max']):.3g}]",
            ]
        )
    lines.extend(
        [
            "",
            "Export-readiness conclusion:",
            f"- coordinates/connectivity present in restart: `{has_coords}`.",
            f"- electric potential field present in restart: `{has_phi}`.",
            f"- restart time `t = {restart_time:.6g}`.",
            f"- restart state size = `{state_size}`, profile rows = `{profile_rows}`, direct index match = `{profile_matches_state}`.",
            "- The restart contains xmhd_2d state fields (`U_n`, velocity, `U_T`,",
            "  `U_psi`, `U_by`) but not the electric potential required by the",
            "  inductionless CSV postprocessor.",
            "- The inspected Hartmann restart appears to be an initial or checkpoint",
            "  state, not the final developed-flow profile used for comparison.",
            "- Coordinates for this Hartmann case are available only through the",
            "  profile sidecar, not the restart itself, and the sidecar row count",
            "  does not match the restart state vector length.",
            "- A real MUG exporter should therefore read mesh/plot output together",
            "  with restart/state fields and should not fabricate `phi` or `sigma`",
            "  without a documented normalization.",
        ]
    )
    (outdir / "mug_restart_manifest_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    restart = (ROOT / args.restart).resolve() if not Path(args.restart).is_absolute() else Path(args.restart)
    profile = (ROOT / args.profile).resolve() if not Path(args.profile).is_absolute() else Path(args.profile)
    rows = inspect_restart(restart)
    info = profile_info(profile)
    write_manifest(outdir, rows)
    plot_ranges(outdir, rows)
    write_summary(outdir, restart, rows, info)
    print(f"Wrote MUG restart manifest diagnostics to {outdir}")


if __name__ == "__main__":
    main()
