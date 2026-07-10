#!/usr/bin/env python3
"""Reduced biased-interface current feedback scan.

Floating sheath loads have zero net current by default. This diagnostic explores
a deliberately biased/electrode-style case in which a prescribed fraction of
the incident ion-saturation current closes through a liquid-metal layer. It
estimates ohmic surface heating and sheet-current magnetic feedback. It is a
triage metric, not a TokaMaker equilibrium response calculation.
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
        classify_feedback_ratio,
        feedback_ratio,
        net_surface_heat_flux,
        plasma_sheath_surface_loads,
        representative_liquid_metal_thermal_properties,
        representative_liquid_metals,
        semi_infinite_constant_flux_temperature_rise,
        sheet_delta_b,
        surface_current_ohmic_heat_flux,
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
    classify_feedback_ratio = coupling.classify_feedback_ratio
    feedback_ratio = coupling.feedback_ratio
    representative_liquid_metals = coupling.representative_liquid_metals
    sheet_delta_b = coupling.sheet_delta_b
    net_surface_heat_flux = plasma_interface.net_surface_heat_flux
    plasma_sheath_surface_loads = plasma_interface.plasma_sheath_surface_loads
    representative_liquid_metal_thermal_properties = plasma_interface.representative_liquid_metal_thermal_properties
    semi_infinite_constant_flux_temperature_rise = plasma_interface.semi_infinite_constant_flux_temperature_rise
    surface_current_ohmic_heat_flux = plasma_interface.surface_current_ohmic_heat_flux


CLASS_VALUE = {
    "weak-feedback": 0,
    "two-way-candidate": 1,
    "strong-two-way-candidate": 2,
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
    parser.add_argument("--bias-fraction", type=float, default=0.02, help="net current as fraction of ion-saturation current")
    parser.add_argument("--conducting-thickness", type=float, default=0.02, help="liquid-metal conducting thickness [m]")
    parser.add_argument("--reference-field", type=float, default=5.0, help="magnetic feedback normalization field [T]")
    parser.add_argument("--heat-transmission", type=float, default=7.0, help="sheath heat transmission coefficient")
    parser.add_argument("--exposure-time", type=float, default=1.0e-2, help="heat-pulse duration [s]")
    return parser.parse_args()


def classify_grid(feedback: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    labels = np.empty(feedback.shape, dtype=object)
    values = np.zeros(feedback.shape, dtype=np.float64)
    for index in np.ndindex(feedback.shape):
        label = classify_feedback_ratio(float(feedback[index]))
        labels[index] = label
        values[index] = CLASS_VALUE[label]
    return labels, values


def compute_scan(args: argparse.Namespace) -> dict[str, np.ndarray | str | float]:
    if args.bias_fraction < 0.0:
        raise ValueError("bias_fraction must be non-negative")
    props = representative_liquid_metals()[args.material]
    thermal_props = representative_liquid_metal_thermal_properties()[args.material]
    densities = np.geomspace(args.density_min, args.density_max, args.n)
    temperatures = np.geomspace(args.temperature_min, args.temperature_max, args.n)
    density_grid, temperature_grid = np.meshgrid(densities, temperatures, indexing="xy")
    loads = plasma_sheath_surface_loads(
        density_grid,
        temperature_grid,
        heat_transmission_coefficient=args.heat_transmission,
    )
    net_current_density = args.bias_fraction * loads.ion_current_density
    ohmic_heat = surface_current_ohmic_heat_flux(
        net_current_density,
        props.electrical_conductivity,
        args.conducting_thickness,
    )
    retained_heat = net_surface_heat_flux(loads.heat_flux + ohmic_heat)
    delta_b = sheet_delta_b(net_current_density, args.conducting_thickness)
    feedback = feedback_ratio(delta_b, args.reference_field)
    class_grid, class_value = classify_grid(np.asarray(feedback))
    semi_infinite_delta_t = semi_infinite_constant_flux_temperature_rise(
        retained_heat,
        args.exposure_time,
        thermal_props.thermal_conductivity,
        thermal_props.density,
        thermal_props.specific_heat,
    )
    return {
        "material": args.material,
        "densities": densities,
        "temperatures": temperatures,
        "density_grid": density_grid,
        "temperature_grid": temperature_grid,
        "ion_current_density": loads.ion_current_density,
        "net_current_density": net_current_density,
        "ohmic_heat_flux": ohmic_heat,
        "plasma_heat_flux": loads.heat_flux,
        "retained_heat_flux": retained_heat,
        "delta_b": delta_b,
        "feedback": feedback,
        "classification": class_grid,
        "classification_value": class_value,
        "semi_infinite_delta_t": semi_infinite_delta_t,
        "bias_fraction": args.bias_fraction,
        "conducting_thickness": args.conducting_thickness,
        "reference_field": args.reference_field,
        "exposure_time": args.exposure_time,
    }


def write_csv(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    with (outdir / "plasma_interface_current_feedback_scan.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "material",
                "electron_density",
                "electron_temperature_ev",
                "ion_current_density",
                "net_current_density",
                "ohmic_heat_flux",
                "plasma_heat_flux",
                "retained_heat_flux",
                "delta_b",
                "feedback",
                "classification",
                "semi_infinite_delta_t",
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
                    "ion_current_density": float(np.asarray(scan["ion_current_density"])[index]),
                    "net_current_density": float(np.asarray(scan["net_current_density"])[index]),
                    "ohmic_heat_flux": float(np.asarray(scan["ohmic_heat_flux"])[index]),
                    "plasma_heat_flux": float(np.asarray(scan["plasma_heat_flux"])[index]),
                    "retained_heat_flux": float(np.asarray(scan["retained_heat_flux"])[index]),
                    "delta_b": float(np.asarray(scan["delta_b"])[index]),
                    "feedback": float(np.asarray(scan["feedback"])[index]),
                    "classification": np.asarray(scan["classification"], dtype=object)[index],
                    "semi_infinite_delta_t": float(np.asarray(scan["semi_infinite_delta_t"])[index]),
                }
            )


def plot_scan(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    panels = [
        ("feedback", "log10 |delta B|/B_ref", "viridis"),
        ("ohmic_heat_flux", "log10 ohmic heat flux [W/m^2]", "magma"),
        ("semi_infinite_delta_t", "log10 surface Delta T [K]", "plasma"),
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
    fig.suptitle(f"{scan['material']} biased-current feedback, bias_fraction={float(scan['bias_fraction']):.3g}")
    fig.tight_layout()
    fig.savefig(outdir / "plasma_interface_current_feedback_scan.png", dpi=180)
    plt.close(fig)


def write_summary(outdir: Path, scan: dict[str, np.ndarray | str | float]) -> None:
    density_target = 1.0e19
    temperature_target = 20.0
    densities = np.asarray(scan["densities"])
    temperatures = np.asarray(scan["temperatures"])
    i = int(np.argmin(np.abs(temperatures - temperature_target)))
    j = int(np.argmin(np.abs(densities - density_target)))
    classification = np.asarray(scan["classification"], dtype=object)
    with (outdir / "plasma_interface_current_feedback_scan_summary.md").open("w", newline="\n") as handle:
        handle.write("# Plasma/LM Biased-Current Feedback Scan\n\n")
        handle.write("This is a reduced biased-interface triage diagnostic, not a TokaMaker response solve.\n\n")
        handle.write("| quantity | selected value |\n")
        handle.write("|---|---:|\n")
        handle.write(f"| material | `{scan['material']}` |\n")
        handle.write(f"| electron density | `{densities[j]:.3e} m^-3` |\n")
        handle.write(f"| electron temperature | `{temperatures[i]:.3e} eV` |\n")
        handle.write(f"| bias fraction | `{float(scan['bias_fraction']):.3e}` |\n")
        handle.write(f"| conducting thickness | `{float(scan['conducting_thickness']):.3e} m` |\n")
        handle.write(f"| net current density | `{np.asarray(scan['net_current_density'])[i, j]:.3e} A/m^2` |\n")
        handle.write(f"| ohmic heat flux | `{np.asarray(scan['ohmic_heat_flux'])[i, j]:.3e} W/m^2` |\n")
        handle.write(f"| feedback ratio | `{np.asarray(scan['feedback'])[i, j]:.3e}` |\n")
        handle.write(f"| feedback classification | `{classification[i, j]}` |\n")
        handle.write(f"| semi-infinite surface Delta T | `{np.asarray(scan['semi_infinite_delta_t'])[i, j]:.3e} K` |\n")
        handle.write("\nFloating sheath cases should use `bias_fraction = 0`; this scan is for biased/electrode-style triage.\n")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scan = compute_scan(args)
    write_csv(outdir, scan)
    plot_scan(outdir, scan)
    write_summary(outdir, scan)
    print(f"Wrote plasma/LM biased-current feedback scan to {outdir}")


if __name__ == "__main__":
    main()
