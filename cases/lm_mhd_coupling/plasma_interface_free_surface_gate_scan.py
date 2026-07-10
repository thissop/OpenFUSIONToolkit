#!/usr/bin/env python3
"""Linear free-surface tractability gate for plasma-loaded LM layers.

This diagnostic combines reduced sheath momentum loading with linear
capillary-gravity-magnetic response scales. It is meant to decide what a future
two-way plasma/liquid-metal calculation would have to resolve; it is not a
nonlinear free-surface simulation.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
try:
    import smplotlib  # noqa: F401
except Exception:  # pragma: no cover
    pass
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "python"))

try:
    from OpenFUSIONToolkit.LM_MHD import (  # noqa: E402
        free_surface_linear_metrics,
        plasma_sheath_surface_loads,
        representative_liquid_metal_thermal_properties,
        representative_liquid_metals,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    interface_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "plasma_interface.py"
    coupling_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    interface_spec = importlib.util.spec_from_file_location("lm_mhd_plasma_interface", interface_path)
    plasma_interface = importlib.util.module_from_spec(interface_spec)
    sys.modules[interface_spec.name] = plasma_interface
    interface_spec.loader.exec_module(plasma_interface)
    coupling_spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(coupling_spec)
    sys.modules[coupling_spec.name] = coupling
    coupling_spec.loader.exec_module(coupling)
    free_surface_linear_metrics = plasma_interface.free_surface_linear_metrics
    plasma_sheath_surface_loads = plasma_interface.plasma_sheath_surface_loads
    representative_liquid_metal_thermal_properties = plasma_interface.representative_liquid_metal_thermal_properties
    representative_liquid_metals = coupling.representative_liquid_metals


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--material", choices=("Li", "PbLi"), default="Li", help="representative liquid metal")
    parser.add_argument("--density-min", type=float, default=1.0e18, help="minimum electron density [m^-3]")
    parser.add_argument("--density-max", type=float, default=1.0e20, help="maximum electron density [m^-3]")
    parser.add_argument("--temperature", type=float, default=20.0, help="electron temperature [eV]")
    parser.add_argument("--wavelength-min", type=float, default=1.0e-3, help="minimum surface wavelength [m]")
    parser.add_argument("--wavelength-max", type=float, default=2.0e-1, help="maximum surface wavelength [m]")
    parser.add_argument("--n", type=int, default=80, help="points per scan axis")
    parser.add_argument("--heat-transmission", type=float, default=7.0, help="sheath heat transmission coefficient")
    parser.add_argument("--layer-thickness", type=float, default=2.0e-3, help="liquid-layer thickness [m]")
    parser.add_argument("--exposure-time", type=float, default=1.0e-2, help="plasma exposure duration [s]")
    parser.add_argument("--magnetic-field", type=float, default=5.0, help="magnetic-field scale [T]")
    parser.add_argument("--magnetic-stiffness-coeff", type=float, default=1.0, help="coefficient multiplying k B^2/mu")
    parser.add_argument("--surface-grid-spacing", type=float, default=1.0e-4, help="candidate surface grid spacing [m]")
    parser.add_argument("--normal-grid-spacing", type=float, default=2.0e-5, help="candidate normal grid spacing [m]")
    parser.add_argument("--min-cells-per-wavelength", type=float, default=32.0, help="surface resolution gate")
    parser.add_argument("--min-cells-across-layer", type=float, default=40.0, help="normal resolution gate")
    parser.add_argument("--linear-displacement-limit", type=float, default=0.05, help="max |eta|/h for linear response")
    parser.add_argument("--timestep-safety", type=float, default=0.05, help="explicit time-step safety factor")
    parser.add_argument("--practical-dt-min", type=float, default=1.0e-7, help="minimum practical explicit step [s]")
    return parser.parse_args()


def compute_scan(args: argparse.Namespace) -> dict[str, np.ndarray | str | float]:
    lm_props = representative_liquid_metals()[args.material]
    thermal_props = representative_liquid_metal_thermal_properties()[args.material]
    densities = np.geomspace(args.density_min, args.density_max, args.n)
    wavelengths = np.geomspace(args.wavelength_min, args.wavelength_max, args.n)
    density_grid, wavelength_grid = np.meshgrid(densities, wavelengths, indexing="xy")
    wavenumber_grid = 2.0 * np.pi / wavelength_grid
    loads = plasma_sheath_surface_loads(
        density_grid,
        args.temperature,
        heat_transmission_coefficient=args.heat_transmission,
    )
    metrics = free_surface_linear_metrics(
        loads.normal_stress,
        wavenumber_grid,
        lm_props.density,
        thermal_props.surface_tension,
        args.layer_thickness,
        args.exposure_time,
        thermal_props.thermal_diffusivity,
        magnetic_field=args.magnetic_field,
        electrical_conductivity=lm_props.electrical_conductivity,
        magnetic_stiffness_coefficient=args.magnetic_stiffness_coeff,
        timestep_safety_factor=args.timestep_safety,
    )
    cells_per_wavelength = wavelength_grid / args.surface_grid_spacing
    cells_across_layer = args.layer_thickness / args.normal_grid_spacing
    linear_ok = metrics.displacement_fraction <= args.linear_displacement_limit
    surface_resolved = cells_per_wavelength >= args.min_cells_per_wavelength
    normal_resolved = cells_across_layer >= args.min_cells_across_layer
    timestep_ok = metrics.explicit_timestep >= args.practical_dt_min
    gate_score = linear_ok.astype(int) + surface_resolved.astype(int) + int(normal_resolved) + timestep_ok.astype(int)
    return {
        "material": args.material,
        "densities": densities,
        "wavelengths": wavelengths,
        "density_grid": density_grid,
        "wavelength_grid": wavelength_grid,
        "wavenumber_grid": wavenumber_grid,
        "normal_stress": loads.normal_stress,
        "heat_flux": loads.heat_flux,
        "surface_stiffness": metrics.stiffness,
        "magnetic_stiffness": metrics.magnetic_stiffness,
        "angular_frequency": metrics.angular_frequency,
        "period": metrics.period,
        "magnetic_damping_rate": metrics.magnetic_damping_rate,
        "magnetic_damping_ratio": metrics.magnetic_damping_ratio,
        "explicit_timestep": metrics.explicit_timestep,
        "static_displacement": metrics.static_displacement,
        "displacement_fraction": metrics.displacement_fraction,
        "thermal_penetration_fraction": metrics.thermal_penetration_fraction,
        "cells_per_wavelength": cells_per_wavelength,
        "cells_across_layer": cells_across_layer,
        "gate_score": gate_score,
        "linear_ok": linear_ok,
        "surface_resolved": surface_resolved,
        "normal_resolved": normal_resolved,
        "timestep_ok": timestep_ok,
        "temperature": args.temperature,
        "layer_thickness": args.layer_thickness,
        "exposure_time": args.exposure_time,
        "magnetic_field": args.magnetic_field,
        "surface_grid_spacing": args.surface_grid_spacing,
        "normal_grid_spacing": args.normal_grid_spacing,
        "min_cells_per_wavelength": args.min_cells_per_wavelength,
        "min_cells_across_layer": args.min_cells_across_layer,
        "linear_displacement_limit": args.linear_displacement_limit,
        "practical_dt_min": args.practical_dt_min,
    }


def write_csv(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    with (outdir / "plasma_interface_free_surface_gate_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "material",
                "electron_density",
                "electron_temperature_ev",
                "wavelength",
                "wavenumber",
                "normal_stress",
                "heat_flux",
                "surface_stiffness",
                "magnetic_stiffness",
                "angular_frequency",
                "period",
                "magnetic_damping_rate",
                "magnetic_damping_ratio",
                "explicit_timestep",
                "static_displacement",
                "displacement_fraction",
                "thermal_penetration_fraction",
                "cells_per_wavelength",
                "cells_across_layer",
                "gate_score",
                "linear_ok",
                "surface_resolved",
                "normal_resolved",
                "timestep_ok",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        density_grid = np.asarray(scan["density_grid"])
        for index in np.ndindex(density_grid.shape):
            writer.writerow(
                {
                    "material": scan["material"],
                    "electron_density": float(density_grid[index]),
                    "electron_temperature_ev": float(scan["temperature"]),
                    "wavelength": float(np.asarray(scan["wavelength_grid"])[index]),
                    "wavenumber": float(np.asarray(scan["wavenumber_grid"])[index]),
                    "normal_stress": float(np.asarray(scan["normal_stress"])[index]),
                    "heat_flux": float(np.asarray(scan["heat_flux"])[index]),
                    "surface_stiffness": float(np.asarray(scan["surface_stiffness"])[index]),
                    "magnetic_stiffness": float(np.asarray(scan["magnetic_stiffness"])[index]),
                    "angular_frequency": float(np.asarray(scan["angular_frequency"])[index]),
                    "period": float(np.asarray(scan["period"])[index]),
                    "magnetic_damping_rate": float(np.asarray(scan["magnetic_damping_rate"])),
                    "magnetic_damping_ratio": float(np.asarray(scan["magnetic_damping_ratio"])[index]),
                    "explicit_timestep": float(np.asarray(scan["explicit_timestep"])[index]),
                    "static_displacement": float(np.asarray(scan["static_displacement"])[index]),
                    "displacement_fraction": float(np.asarray(scan["displacement_fraction"])[index]),
                    "thermal_penetration_fraction": float(np.asarray(scan["thermal_penetration_fraction"])),
                    "cells_per_wavelength": float(np.asarray(scan["cells_per_wavelength"])[index]),
                    "cells_across_layer": float(scan["cells_across_layer"]),
                    "gate_score": int(np.asarray(scan["gate_score"])[index]),
                    "linear_ok": bool(np.asarray(scan["linear_ok"])[index]),
                    "surface_resolved": bool(np.asarray(scan["surface_resolved"])[index]),
                    "normal_resolved": bool(scan["normal_resolved"]),
                    "timestep_ok": bool(np.asarray(scan["timestep_ok"])[index]),
                }
            )


def plot_scan(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    panels = [
        ("normal_stress", "log10 sheath normal stress [Pa]", "viridis"),
        ("displacement_fraction", "log10 |eta| / h", "magma"),
        ("explicit_timestep", "log10 explicit dt gate [s]", "plasma"),
        ("magnetic_damping_ratio", "log10 magnetic damping ratio", "cividis"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.4), sharex=True, sharey=True)
    densities = np.asarray(scan["densities"])
    wavelengths = np.asarray(scan["wavelengths"])
    for ax, (key, title, cmap) in zip(axes.ravel(), panels):
        data = np.asarray(scan[key])
        image = ax.pcolormesh(
            densities,
            wavelengths,
            np.log10(np.maximum(np.abs(data), np.finfo(float).tiny)),
            shading="auto",
            cmap=cmap,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("electron density [m^-3]")
        ax.set_ylabel("surface wavelength [m]")
        ax.set_title(title)
        fig.colorbar(image, ax=ax)
    fig.suptitle(f"{scan['material']} linear free-surface gate, Te={float(scan['temperature']):.1f} eV")
    fig.tight_layout()
    fig.savefig(outdir / "plasma_interface_free_surface_gate_scan.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    density_target = 1.0e19
    wavelength_target = 2.0e-2
    densities = np.asarray(scan["densities"])
    wavelengths = np.asarray(scan["wavelengths"])
    j = int(np.argmin(np.abs(densities - density_target)))
    i = int(np.argmin(np.abs(wavelengths - wavelength_target)))
    gate_score = np.asarray(scan["gate_score"])
    full_gate_fraction = float(np.count_nonzero(gate_score == 4) / gate_score.size)
    with (outdir / "plasma_interface_free_surface_gate_scan_summary.md").open("w", newline="\n") as handle:
        handle.write("# Plasma/LM Free-Surface Gate Scan\n\n")
        handle.write(
            "This is a linear capillary-gravity-magnetic tractability diagnostic, "
            "not a nonlinear free-surface simulation.\n\n"
        )
        handle.write("| quantity | selected value |\n")
        handle.write("|---|---:|\n")
        handle.write(f"| material | `{scan['material']}` |\n")
        handle.write(f"| electron density | `{densities[j]:.3e} m^-3` |\n")
        handle.write(f"| electron temperature | `{float(scan['temperature']):.3e} eV` |\n")
        handle.write(f"| wavelength | `{wavelengths[i]:.3e} m` |\n")
        handle.write(f"| magnetic field | `{float(scan['magnetic_field']):.3e} T` |\n")
        handle.write(f"| layer thickness | `{float(scan['layer_thickness']):.3e} m` |\n")
        handle.write(f"| normal stress | `{np.asarray(scan['normal_stress'])[i, j]:.3e} Pa` |\n")
        handle.write(f"| surface stiffness | `{np.asarray(scan['surface_stiffness'])[i, j]:.3e} Pa/m` |\n")
        handle.write(f"| magnetic stiffness | `{np.asarray(scan['magnetic_stiffness'])[i, j]:.3e} Pa/m` |\n")
        handle.write(f"| static displacement | `{np.asarray(scan['static_displacement'])[i, j]:.3e} m` |\n")
        handle.write(f"| displacement / layer | `{np.asarray(scan['displacement_fraction'])[i, j]:.3e}` |\n")
        handle.write(f"| wave period | `{np.asarray(scan['period'])[i, j]:.3e} s` |\n")
        handle.write(f"| explicit dt gate | `{np.asarray(scan['explicit_timestep'])[i, j]:.3e} s` |\n")
        handle.write(f"| magnetic damping ratio | `{np.asarray(scan['magnetic_damping_ratio'])[i, j]:.3e}` |\n")
        handle.write(f"| thermal penetration / layer | `{float(np.asarray(scan['thermal_penetration_fraction'])):.3e}` |\n")
        handle.write(f"| cells per wavelength | `{np.asarray(scan['cells_per_wavelength'])[i, j]:.3e}` |\n")
        handle.write(f"| cells across layer | `{float(scan['cells_across_layer']):.3e}` |\n")
        handle.write(f"| full gate pass fraction | `{full_gate_fraction:.3f}` |\n")
        handle.write("\nGate interpretation:\n")
        handle.write("- `linear_ok` requires `|eta|/h` below the configured small-displacement limit.\n")
        handle.write("- `surface_resolved` requires enough cells per surface wavelength.\n")
        handle.write("- `normal_resolved` requires enough cells across the liquid layer.\n")
        handle.write("- `timestep_ok` requires the explicit wave/damping step to exceed the configured practical floor.\n")
        handle.write(
            "\nThe magnetic stiffness scale uses `coefficient * k B^2 / mu`; the coefficient is a model input "
            "because exact prefactors depend on geometry and electromagnetic boundary conditions.\n"
        )
        handle.write("\nVerification backing: `src/tests/physics/test_lm_mhd_plasma_interface.py`.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scan = compute_scan(args)
    write_csv(outdir, scan)
    plot_scan(outdir, scan)
    write_summary(outdir, scan)
    print(f"Wrote plasma/LM free-surface gate scan to {outdir}")


if __name__ == "__main__":
    main()
