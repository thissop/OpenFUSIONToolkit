#!/usr/bin/env python3
"""Reduced plasma/liquid-metal interface sheath-load scan.

This diagnostic maps edge-plasma density and electron temperature to sheath
heat flux, incident ion current, ion momentum load, and a linear static
capillary-gravity displacement scale. It is a boundary-load/regime diagnostic,
not a plasma sheath solver and not a nonlinear free-surface simulation.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
try:
    import smplotlib  # noqa: F401  (applies the style on import)
except Exception:  # pragma: no cover - optional plotting style
    pass
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "python"))

try:
    from OpenFUSIONToolkit.LM_MHD import (  # noqa: E402
        plasma_sheath_surface_loads,
        static_capillary_gravity_displacement,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    interface_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "plasma_interface.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_plasma_interface", interface_path)
    plasma_interface = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = plasma_interface
    spec.loader.exec_module(plasma_interface)
    plasma_sheath_surface_loads = plasma_interface.plasma_sheath_surface_loads
    static_capillary_gravity_displacement = plasma_interface.static_capillary_gravity_displacement


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--density-min", type=float, default=1.0e18, help="minimum electron density [m^-3]")
    parser.add_argument("--density-max", type=float, default=1.0e20, help="maximum electron density [m^-3]")
    parser.add_argument("--temperature-min", type=float, default=2.0, help="minimum electron temperature [eV]")
    parser.add_argument("--temperature-max", type=float, default=100.0, help="maximum electron temperature [eV]")
    parser.add_argument("--n", type=int, default=80, help="points per scan axis")
    parser.add_argument("--heat-transmission", type=float, default=7.0, help="sheath heat transmission coefficient")
    parser.add_argument("--lm-density", type=float, default=500.0, help="liquid-metal density for surface response [kg/m^3]")
    parser.add_argument("--surface-tension", type=float, default=0.38, help="surface tension [N/m]")
    parser.add_argument("--wavelength", type=float, default=0.02, help="surface-response wavelength [m]")
    return parser.parse_args()


def compute_scan(args: argparse.Namespace) -> dict[str, np.ndarray]:
    densities = np.geomspace(args.density_min, args.density_max, args.n)
    temperatures = np.geomspace(args.temperature_min, args.temperature_max, args.n)
    density_grid, temperature_grid = np.meshgrid(densities, temperatures, indexing="xy")
    loads = plasma_sheath_surface_loads(
        density_grid,
        temperature_grid,
        heat_transmission_coefficient=args.heat_transmission,
    )
    wavenumber = 2.0 * np.pi / args.wavelength
    displacement = static_capillary_gravity_displacement(
        loads.normal_stress,
        wavenumber,
        density=args.lm_density,
        surface_tension=args.surface_tension,
    )
    return {
        "densities": densities,
        "temperatures": temperatures,
        "density_grid": density_grid,
        "temperature_grid": temperature_grid,
        "heat_flux": loads.heat_flux,
        "ion_current_density": loads.ion_current_density,
        "normal_stress": loads.normal_stress,
        "surface_displacement": displacement,
    }


def write_csv(outdir: Path, scan: dict[str, np.ndarray]) -> None:
    with (outdir / "plasma_interface_sheath_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "electron_density",
                "electron_temperature_ev",
                "heat_flux",
                "ion_current_density",
                "normal_stress",
                "surface_displacement",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for index in np.ndindex(scan["density_grid"].shape):
            writer.writerow(
                {
                    "electron_density": float(scan["density_grid"][index]),
                    "electron_temperature_ev": float(scan["temperature_grid"][index]),
                    "heat_flux": float(scan["heat_flux"][index]),
                    "ion_current_density": float(scan["ion_current_density"][index]),
                    "normal_stress": float(scan["normal_stress"][index]),
                    "surface_displacement": float(scan["surface_displacement"][index]),
                }
            )


def plot_scan(outdir: Path, scan: dict[str, np.ndarray]) -> None:
    panels = [
        ("heat_flux", "log10 heat flux [W/m^2]"),
        ("ion_current_density", "log10 incident ion current [A/m^2]"),
        ("surface_displacement", "log10 linear displacement [m]"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.2), sharex=True, sharey=True)
    for ax, (key, title) in zip(axes, panels):
        data = scan[key]
        image = ax.pcolormesh(
            scan["densities"],
            scan["temperatures"],
            np.log10(np.maximum(np.abs(data), np.finfo(float).tiny)),
            shading="auto",
            cmap="viridis",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("electron density [m^-3]")
        ax.set_ylabel("electron temperature [eV]")
        ax.set_title(title)
        fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(outdir / "plasma_interface_sheath_scan.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, args: argparse.Namespace, scan: dict[str, np.ndarray]) -> None:
    density_target = 1.0e19
    temperature_target = 20.0
    i = int(np.argmin(np.abs(scan["temperatures"] - temperature_target)))
    j = int(np.argmin(np.abs(scan["densities"] - density_target)))
    with (outdir / "plasma_interface_sheath_scan_summary.md").open("w", newline="\n") as handle:
        handle.write("# Plasma/LM Interface Sheath-Load Scan\n\n")
        handle.write("This is a reduced boundary-load diagnostic, not a sheath solver or free-surface simulation.\n\n")
        handle.write("| quantity | selected value |\n")
        handle.write("|---|---:|\n")
        handle.write(f"| electron density | `{scan['densities'][j]:.3e} m^-3` |\n")
        handle.write(f"| electron temperature | `{scan['temperatures'][i]:.3e} eV` |\n")
        handle.write(f"| heat transmission coefficient | `{args.heat_transmission:.3e}` |\n")
        handle.write(f"| surface-response wavelength | `{args.wavelength:.3e} m` |\n")
        handle.write(f"| heat flux | `{scan['heat_flux'][i, j]:.3e} W/m^2` |\n")
        handle.write(f"| incident ion current density | `{scan['ion_current_density'][i, j]:.3e} A/m^2` |\n")
        handle.write(f"| ion normal stress | `{scan['normal_stress'][i, j]:.3e} Pa` |\n")
        handle.write(f"| linear static displacement scale | `{scan['surface_displacement'][i, j]:.3e} m` |\n")
        handle.write("\nVerification backing: `src/tests/physics/test_lm_mhd_plasma_interface.py`.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scan = compute_scan(args)
    write_csv(outdir, scan)
    plot_scan(outdir, scan)
    write_summary(outdir, args, scan)
    print(f"Wrote plasma/LM sheath-load scan to {outdir}")


if __name__ == "__main__":
    main()
