#!/usr/bin/env python3
"""Reduced thermocapillary surface-flow scan for plasma-heated LM layers.

This diagnostic maps sheath-driven analytic surface temperature rise to a
Marangoni shear stress and the corresponding thin-layer Couette velocity. It
is a linear fallback benchmark, not a nonlinear free-surface calculation.
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
        marangoni_shear_stress,
        net_surface_heat_flux,
        plasma_sheath_surface_loads,
        representative_liquid_metal_thermal_properties,
        representative_liquid_metals,
        semi_infinite_constant_flux_temperature_rise,
        thermocapillary_couette_mean_velocity,
        thermocapillary_couette_profile,
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
    marangoni_shear_stress = plasma_interface.marangoni_shear_stress
    net_surface_heat_flux = plasma_interface.net_surface_heat_flux
    plasma_sheath_surface_loads = plasma_interface.plasma_sheath_surface_loads
    representative_liquid_metal_thermal_properties = plasma_interface.representative_liquid_metal_thermal_properties
    semi_infinite_constant_flux_temperature_rise = plasma_interface.semi_infinite_constant_flux_temperature_rise
    thermocapillary_couette_mean_velocity = plasma_interface.thermocapillary_couette_mean_velocity
    thermocapillary_couette_profile = plasma_interface.thermocapillary_couette_profile
    representative_liquid_metals = coupling.representative_liquid_metals


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="cases/lm_mhd_coupling/results", help="output directory")
    parser.add_argument("--material", choices=("Li", "PbLi"), default="Li", help="representative liquid metal")
    parser.add_argument("--density-min", type=float, default=1.0e18, help="minimum electron density [m^-3]")
    parser.add_argument("--density-max", type=float, default=1.0e20, help="maximum electron density [m^-3]")
    parser.add_argument("--temperature-min", type=float, default=2.0, help="minimum electron temperature [eV]")
    parser.add_argument("--temperature-max", type=float, default=100.0, help="maximum electron temperature [eV]")
    parser.add_argument("--n", type=int, default=80, help="points per scan axis")
    parser.add_argument("--heat-transmission", type=float, default=7.0, help="sheath heat transmission coefficient")
    parser.add_argument("--exposure-time", type=float, default=1.0e-2, help="heat-pulse duration [s]")
    parser.add_argument("--layer-thickness", type=float, default=2.0e-3, help="liquid-layer thickness [m]")
    parser.add_argument("--gradient-length", type=float, default=2.0e-2, help="surface temperature-gradient length [m]")
    parser.add_argument("--d-sigma-dt", type=float, default=-1.0e-4, help="surface-tension temperature coefficient [N/m/K]")
    return parser.parse_args()


def compute_scan(args: argparse.Namespace) -> dict[str, np.ndarray | str | float]:
    lm_props = representative_liquid_metals()[args.material]
    thermal_props = representative_liquid_metal_thermal_properties()[args.material]
    dynamic_viscosity = lm_props.density * lm_props.kinematic_viscosity
    densities = np.geomspace(args.density_min, args.density_max, args.n)
    temperatures = np.geomspace(args.temperature_min, args.temperature_max, args.n)
    density_grid, temperature_grid = np.meshgrid(densities, temperatures, indexing="xy")
    loads = plasma_sheath_surface_loads(
        density_grid,
        temperature_grid,
        heat_transmission_coefficient=args.heat_transmission,
    )
    retained_heat = net_surface_heat_flux(loads.heat_flux)
    surface_delta_t = semi_infinite_constant_flux_temperature_rise(
        retained_heat,
        args.exposure_time,
        thermal_props.thermal_conductivity,
        thermal_props.density,
        thermal_props.specific_heat,
    )
    surface_temperature_gradient = surface_delta_t / args.gradient_length
    shear = marangoni_shear_stress(surface_temperature_gradient, args.d_sigma_dt)
    surface_velocity = thermocapillary_couette_profile(
        args.layer_thickness,
        args.layer_thickness,
        dynamic_viscosity,
        surface_temperature_gradient,
        args.d_sigma_dt,
    )
    mean_velocity = thermocapillary_couette_mean_velocity(
        args.layer_thickness,
        dynamic_viscosity,
        surface_temperature_gradient,
        args.d_sigma_dt,
    )
    return {
        "material": args.material,
        "densities": densities,
        "temperatures": temperatures,
        "density_grid": density_grid,
        "temperature_grid": temperature_grid,
        "heat_flux": loads.heat_flux,
        "surface_delta_t": surface_delta_t,
        "surface_temperature_gradient": surface_temperature_gradient,
        "marangoni_shear": shear,
        "surface_velocity": surface_velocity,
        "mean_velocity": mean_velocity,
        "layer_thickness": args.layer_thickness,
        "gradient_length": args.gradient_length,
        "d_sigma_dt": args.d_sigma_dt,
        "dynamic_viscosity": dynamic_viscosity,
        "exposure_time": args.exposure_time,
    }


def write_csv(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    with (outdir / "plasma_interface_marangoni_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "material",
                "electron_density",
                "electron_temperature_ev",
                "heat_flux",
                "surface_delta_t",
                "surface_temperature_gradient",
                "marangoni_shear",
                "surface_velocity",
                "mean_velocity",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        density_grid = np.asarray(scan["density_grid"])
        temperature_grid = np.asarray(scan["temperature_grid"])
        for index in np.ndindex(density_grid.shape):
            writer.writerow(
                {
                    "material": scan["material"],
                    "electron_density": float(density_grid[index]),
                    "electron_temperature_ev": float(temperature_grid[index]),
                    "heat_flux": float(np.asarray(scan["heat_flux"])[index]),
                    "surface_delta_t": float(np.asarray(scan["surface_delta_t"])[index]),
                    "surface_temperature_gradient": float(np.asarray(scan["surface_temperature_gradient"])[index]),
                    "marangoni_shear": float(np.asarray(scan["marangoni_shear"])[index]),
                    "surface_velocity": float(np.asarray(scan["surface_velocity"])[index]),
                    "mean_velocity": float(np.asarray(scan["mean_velocity"])[index]),
                }
            )


def plot_scan(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    panels = [
        ("surface_delta_t", "log10 surface Delta T [K]", "magma"),
        ("marangoni_shear", "log10 |Marangoni shear| [Pa]", "viridis"),
        ("surface_velocity", "log10 |surface velocity| [m/s]", "plasma"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.2), sharex=True, sharey=True)
    densities = np.asarray(scan["densities"])
    temperatures = np.asarray(scan["temperatures"])
    for ax, (key, title, cmap) in zip(axes, panels):
        data = np.asarray(scan[key])
        image = ax.pcolormesh(
            densities,
            temperatures,
            np.log10(np.maximum(np.abs(data), np.finfo(float).tiny)),
            shading="auto",
            cmap=cmap,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("electron density [m^-3]")
        ax.set_ylabel("electron temperature [eV]")
        ax.set_title(title)
        fig.colorbar(image, ax=ax)
    fig.suptitle(f"{scan['material']} thermocapillary fallback scan")
    fig.tight_layout()
    fig.savefig(outdir / "plasma_interface_marangoni_scan.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    density_target = 1.0e19
    temperature_target = 20.0
    densities = np.asarray(scan["densities"])
    temperatures = np.asarray(scan["temperatures"])
    i = int(np.argmin(np.abs(temperatures - temperature_target)))
    j = int(np.argmin(np.abs(densities - density_target)))
    with (outdir / "plasma_interface_marangoni_scan_summary.md").open("w", newline="\n") as handle:
        handle.write("# Plasma/LM Thermocapillary Fallback Scan\n\n")
        handle.write("This is a linear Marangoni/Couette diagnostic, not a nonlinear free-surface simulation.\n\n")
        handle.write("| quantity | selected value |\n")
        handle.write("|---|---:|\n")
        handle.write(f"| material | `{scan['material']}` |\n")
        handle.write(f"| electron density | `{densities[j]:.3e} m^-3` |\n")
        handle.write(f"| electron temperature | `{temperatures[i]:.3e} eV` |\n")
        handle.write(f"| layer thickness | `{float(scan['layer_thickness']):.3e} m` |\n")
        handle.write(f"| gradient length | `{float(scan['gradient_length']):.3e} m` |\n")
        handle.write(f"| d sigma / dT | `{float(scan['d_sigma_dt']):.3e} N/m/K` |\n")
        handle.write(f"| surface Delta T | `{np.asarray(scan['surface_delta_t'])[i, j]:.3e} K` |\n")
        handle.write(f"| Marangoni shear | `{np.asarray(scan['marangoni_shear'])[i, j]:.3e} Pa` |\n")
        handle.write(f"| surface velocity | `{np.asarray(scan['surface_velocity'])[i, j]:.3e} m/s` |\n")
        handle.write(f"| mean velocity | `{np.asarray(scan['mean_velocity'])[i, j]:.3e} m/s` |\n")
        handle.write("\nVerification backing: `src/tests/physics/test_lm_mhd_plasma_interface.py`.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scan = compute_scan(args)
    write_csv(outdir, scan)
    plot_scan(outdir, scan)
    write_summary(outdir, scan)
    print(f"Wrote plasma/LM thermocapillary scan to {outdir}")


if __name__ == "__main__":
    main()
