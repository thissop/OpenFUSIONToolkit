#!/usr/bin/env python3
"""Reduced plasma/liquid-metal surface thermal-response scan.

This diagnostic composes the sheath heat-flux contract with two analytic
thermal responses: a well-mixed surface layer and a semi-infinite
constant-flux conduction limit. It is a reference/triage calculation, not a
thermal CFD or free-surface simulation.
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
        classify_thermal_response_regime,
        lumped_surface_temperature_rise,
        net_surface_heat_flux,
        plasma_sheath_surface_loads,
        representative_liquid_metal_thermal_properties,
        semi_infinite_constant_flux_temperature_rise,
    )
except FileNotFoundError:  # pragma: no cover - source-tree convenience path
    import importlib.util  # noqa: E402

    interface_path = ROOT / "src" / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "plasma_interface.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_plasma_interface", interface_path)
    plasma_interface = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = plasma_interface
    spec.loader.exec_module(plasma_interface)
    classify_thermal_response_regime = plasma_interface.classify_thermal_response_regime
    lumped_surface_temperature_rise = plasma_interface.lumped_surface_temperature_rise
    net_surface_heat_flux = plasma_interface.net_surface_heat_flux
    plasma_sheath_surface_loads = plasma_interface.plasma_sheath_surface_loads
    representative_liquid_metal_thermal_properties = plasma_interface.representative_liquid_metal_thermal_properties
    semi_infinite_constant_flux_temperature_rise = plasma_interface.semi_infinite_constant_flux_temperature_rise


REGIME_VALUE = {
    "semi-infinite": 0,
    "transient-slab": 1,
    "through-thickness": 2,
}


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
    parser.add_argument("--layer-thickness", type=float, default=2.0e-3, help="well-mixed surface-layer thickness [m]")
    parser.add_argument("--evaporation-mass-flux", type=float, default=0.0, help="constant evaporative mass flux [kg/m^2/s]")
    return parser.parse_args()


def compute_scan(args: argparse.Namespace) -> dict[str, np.ndarray | str | float]:
    props = representative_liquid_metal_thermal_properties()[args.material]
    densities = np.geomspace(args.density_min, args.density_max, args.n)
    temperatures = np.geomspace(args.temperature_min, args.temperature_max, args.n)
    density_grid, temperature_grid = np.meshgrid(densities, temperatures, indexing="xy")
    loads = plasma_sheath_surface_loads(
        density_grid,
        temperature_grid,
        heat_transmission_coefficient=args.heat_transmission,
    )
    latent_heat = 0.0 if props.latent_heat is None else props.latent_heat
    retained_heat = net_surface_heat_flux(
        loads.heat_flux,
        evaporation_mass_flux=args.evaporation_mass_flux,
        latent_heat=latent_heat,
    )
    semi_infinite_delta_t = semi_infinite_constant_flux_temperature_rise(
        retained_heat,
        args.exposure_time,
        props.thermal_conductivity,
        props.density,
        props.specific_heat,
    )
    lumped_delta_t = lumped_surface_temperature_rise(
        retained_heat,
        args.exposure_time,
        args.layer_thickness,
        props.density,
        props.specific_heat,
    )
    regime = classify_thermal_response_regime(
        args.exposure_time,
        args.layer_thickness,
        props.thermal_conductivity,
        props.density,
        props.specific_heat,
    )
    return {
        "material": args.material,
        "densities": densities,
        "temperatures": temperatures,
        "density_grid": density_grid,
        "temperature_grid": temperature_grid,
        "heat_flux": loads.heat_flux,
        "retained_heat_flux": retained_heat,
        "semi_infinite_delta_t": semi_infinite_delta_t,
        "lumped_delta_t": lumped_delta_t,
        "thermal_regime": str(regime),
        "exposure_time": args.exposure_time,
        "layer_thickness": args.layer_thickness,
        "density": props.density,
        "specific_heat": props.specific_heat,
        "thermal_conductivity": props.thermal_conductivity,
    }


def write_csv(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    with (outdir / "plasma_interface_thermal_response_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "material",
                "electron_density",
                "electron_temperature_ev",
                "heat_flux",
                "retained_heat_flux",
                "semi_infinite_delta_t",
                "lumped_delta_t",
                "thermal_regime",
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
                    "retained_heat_flux": float(np.asarray(scan["retained_heat_flux"])[index]),
                    "semi_infinite_delta_t": float(np.asarray(scan["semi_infinite_delta_t"])[index]),
                    "lumped_delta_t": float(np.asarray(scan["lumped_delta_t"])[index]),
                    "thermal_regime": scan["thermal_regime"],
                }
            )


def plot_scan(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    panels = [
        ("retained_heat_flux", "log10 retained heat flux [W/m^2]"),
        ("semi_infinite_delta_t", "log10 semi-infinite Delta T [K]"),
        ("lumped_delta_t", "log10 lumped-layer Delta T [K]"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.2), sharex=True, sharey=True)
    densities = np.asarray(scan["densities"])
    temperatures = np.asarray(scan["temperatures"])
    for ax, (key, title) in zip(axes, panels):
        data = np.asarray(scan[key])
        image = ax.pcolormesh(
            densities,
            temperatures,
            np.log10(np.maximum(np.abs(data), np.finfo(float).tiny)),
            shading="auto",
            cmap="magma",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("electron density [m^-3]")
        ax.set_ylabel("electron temperature [eV]")
        ax.set_title(title)
        fig.colorbar(image, ax=ax)
    fig.suptitle(f"{scan['material']} thermal response, regime: {scan['thermal_regime']}")
    fig.tight_layout()
    fig.savefig(outdir / "plasma_interface_thermal_response_scan.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    density_target = 1.0e19
    temperature_target = 20.0
    densities = np.asarray(scan["densities"])
    temperatures = np.asarray(scan["temperatures"])
    i = int(np.argmin(np.abs(temperatures - temperature_target)))
    j = int(np.argmin(np.abs(densities - density_target)))
    with (outdir / "plasma_interface_thermal_response_scan_summary.md").open("w", newline="\n") as handle:
        handle.write("# Plasma/LM Thermal-Response Scan\n\n")
        handle.write("This is a reduced analytic thermal-response diagnostic, not a thermal CFD result.\n\n")
        handle.write("| quantity | selected value |\n")
        handle.write("|---|---:|\n")
        handle.write(f"| material | `{scan['material']}` |\n")
        handle.write(f"| electron density | `{densities[j]:.3e} m^-3` |\n")
        handle.write(f"| electron temperature | `{temperatures[i]:.3e} eV` |\n")
        handle.write(f"| exposure time | `{float(scan['exposure_time']):.3e} s` |\n")
        handle.write(f"| layer thickness | `{float(scan['layer_thickness']):.3e} m` |\n")
        handle.write(f"| thermal regime | `{scan['thermal_regime']}` |\n")
        handle.write(f"| heat flux | `{np.asarray(scan['heat_flux'])[i, j]:.3e} W/m^2` |\n")
        handle.write(f"| retained heat flux | `{np.asarray(scan['retained_heat_flux'])[i, j]:.3e} W/m^2` |\n")
        handle.write(f"| semi-infinite surface Delta T | `{np.asarray(scan['semi_infinite_delta_t'])[i, j]:.3e} K` |\n")
        handle.write(f"| lumped-layer Delta T | `{np.asarray(scan['lumped_delta_t'])[i, j]:.3e} K` |\n")
        handle.write("\nVerification backing: `src/tests/physics/test_lm_mhd_plasma_interface.py`.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scan = compute_scan(args)
    write_csv(outdir, scan)
    plot_scan(outdir, scan)
    write_summary(outdir, scan)
    print(f"Wrote plasma/LM thermal-response scan to {outdir}")


if __name__ == "__main__":
    main()
