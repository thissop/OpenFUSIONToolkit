#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Reference utilities for plasma/liquid-metal interface load accounting.

These helpers provide formula-level boundary-load contracts for future
plasma/liquid-metal coupling work. They are not a plasma sheath solver and do
not advance a liquid free surface. The sign convention is explicit: positive
heat flux, particle flux, current density, and normal stress are loads delivered
to the liquid-metal surface along the chosen surface normal.
'''
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

E_CHARGE = 1.602176634e-19
ELECTRON_MASS = 9.1093837015e-31
PROTON_MASS = 1.67262192369e-27
DEUTERIUM_ION_MASS = 2.013553212745 * PROTON_MASS
STANDARD_GRAVITY = 9.80665
VACUUM_PERMEABILITY = 4.0e-7 * np.pi


@dataclass(frozen=True)
class PlasmaLiquidInterfaceLoads:
    '''Surface loads from a reduced plasma/sheath model.

    All fields use SI units and are positive into the liquid-metal surface under
    the caller's chosen normal convention. `ion_current_density` is the incident
    ion-current magnitude. `net_current_density` is the electrical current
    actually delivered to the liquid metal; it defaults to zero for a floating
    interface and should be supplied explicitly for biased/electrode cases.
    '''

    particle_flux: np.ndarray
    incident_mass_flux: np.ndarray
    ion_current_density: np.ndarray
    net_current_density: np.ndarray
    heat_flux: np.ndarray
    normal_stress: np.ndarray


@dataclass(frozen=True)
class LiquidMetalThermalProperties:
    '''Uniform thermal properties for reduced surface-response estimates.'''

    name: str
    density: float
    specific_heat: float
    thermal_conductivity: float
    surface_tension: float | None = None
    latent_heat: float | None = None
    note: str = ""

    def __post_init__(self) -> None:
        _require_positive("density", self.density)
        _require_positive("specific_heat", self.specific_heat)
        _require_positive("thermal_conductivity", self.thermal_conductivity)
        if self.surface_tension is not None:
            _require_positive("surface_tension", self.surface_tension)
        if self.latent_heat is not None:
            _require_positive("latent_heat", self.latent_heat)

    @property
    def thermal_diffusivity(self) -> float:
        '''Return alpha = k / (rho cp).'''

        return thermal_diffusivity(self.thermal_conductivity, self.density, self.specific_heat)


@dataclass(frozen=True)
class FreeSurfaceLinearMetrics:
    '''Linear free-surface response scales for a plasma-loaded LM interface.

    These are pre-solve gate metrics for deciding whether a moving-surface
    calculation is likely to be stiff, underresolved, or outside the small-slope
    limit. They are not a nonlinear free-surface model.
    '''

    wavenumber: np.ndarray
    wavelength: np.ndarray
    stiffness: np.ndarray
    magnetic_stiffness: np.ndarray
    angular_frequency: np.ndarray
    period: np.ndarray
    magnetic_damping_rate: np.ndarray
    magnetic_damping_ratio: np.ndarray
    explicit_timestep: np.ndarray
    static_displacement: np.ndarray
    displacement_fraction: np.ndarray
    thermal_penetration_fraction: np.ndarray


def representative_liquid_metal_thermal_properties() -> dict[str, LiquidMetalThermalProperties]:
    '''Return rounded Li and PbLi thermal scales for exploratory diagnostics.

    These values are intentionally approximate. Temperature-dependent property
    correlations should replace them before any design or predictive claim.
    '''

    return {
        "Li": LiquidMetalThermalProperties(
            name="Li, representative near 500 C",
            density=5.0e2,
            specific_heat=4.2e3,
            thermal_conductivity=5.0e1,
            surface_tension=0.38,
            latent_heat=2.0e7,
            note="Rounded liquid-lithium scale for verification/regime scans.",
        ),
        "PbLi": LiquidMetalThermalProperties(
            name="PbLi, representative near 500 C",
            density=9.5e3,
            specific_heat=1.5e2,
            thermal_conductivity=1.6e1,
            surface_tension=0.45,
            latent_heat=8.0e5,
            note="Rounded PbLi scale for verification/regime scans.",
        ),
    }


def bohm_sound_speed(
    electron_temperature_ev: float | np.ndarray,
    ion_mass: float = DEUTERIUM_ION_MASS,
    ion_temperature_ev: float | np.ndarray = 0.0,
    ion_gamma: float = 1.0,
    ion_charge_state: float = 1.0,
) -> float | np.ndarray:
    '''Return the Bohm sound speed c_s = sqrt(e*(Z*T_e + gamma_i*T_i)/m_i).'''

    _require_positive("ion_mass", ion_mass)
    _require_nonnegative("ion_gamma", ion_gamma)
    _require_positive("ion_charge_state", ion_charge_state)
    te = _positive_array("electron_temperature_ev", electron_temperature_ev)
    ti = _nonnegative_array("ion_temperature_ev", ion_temperature_ev)
    speed = np.sqrt(E_CHARGE * (ion_charge_state * te + ion_gamma * ti) / ion_mass)
    return _maybe_scalar(speed)


def sheath_particle_flux(
    electron_density: float | np.ndarray,
    electron_temperature_ev: float | np.ndarray,
    ion_mass: float = DEUTERIUM_ION_MASS,
    ion_temperature_ev: float | np.ndarray = 0.0,
    ion_gamma: float = 1.0,
    ion_charge_state: float = 1.0,
    sheath_density_factor: float = 1.0,
) -> float | np.ndarray:
    '''Return incident ion particle flux Gamma_i = f_s n_e c_s.'''

    _require_nonnegative("sheath_density_factor", sheath_density_factor)
    density = _positive_array("electron_density", electron_density)
    sound_speed = bohm_sound_speed(
        electron_temperature_ev,
        ion_mass=ion_mass,
        ion_temperature_ev=ion_temperature_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=ion_charge_state,
    )
    return _maybe_scalar(sheath_density_factor * density * np.asarray(sound_speed))


def ion_saturation_current_density(
    electron_density: float | np.ndarray,
    electron_temperature_ev: float | np.ndarray,
    ion_mass: float = DEUTERIUM_ION_MASS,
    ion_temperature_ev: float | np.ndarray = 0.0,
    ion_gamma: float = 1.0,
    ion_charge_state: float = 1.0,
    sheath_density_factor: float = 1.0,
) -> float | np.ndarray:
    '''Return incident ion-current density j_i = Z e Gamma_i.'''

    gamma_i = sheath_particle_flux(
        electron_density,
        electron_temperature_ev,
        ion_mass=ion_mass,
        ion_temperature_ev=ion_temperature_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=ion_charge_state,
        sheath_density_factor=sheath_density_factor,
    )
    return _maybe_scalar(ion_charge_state * E_CHARGE * np.asarray(gamma_i))


def sheath_heat_transmission_flux(
    electron_density: float | np.ndarray,
    electron_temperature_ev: float | np.ndarray,
    ion_mass: float = DEUTERIUM_ION_MASS,
    ion_temperature_ev: float | np.ndarray = 0.0,
    ion_gamma: float = 1.0,
    ion_charge_state: float = 1.0,
    sheath_density_factor: float = 1.0,
    heat_transmission_coefficient: float = 7.0,
) -> float | np.ndarray:
    '''Return q_parallel = gamma_sh Gamma_i e T_e.

    The coefficient is model input, not a universal constant. A value near 7 is
    a common deuterium edge-plasma scale and is useful for first-pass load
    accounting, but production cases should supply the sheath model's value.
    '''

    _require_positive("heat_transmission_coefficient", heat_transmission_coefficient)
    gamma_i = sheath_particle_flux(
        electron_density,
        electron_temperature_ev,
        ion_mass=ion_mass,
        ion_temperature_ev=ion_temperature_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=ion_charge_state,
        sheath_density_factor=sheath_density_factor,
    )
    te = _positive_array("electron_temperature_ev", electron_temperature_ev)
    return _maybe_scalar(heat_transmission_coefficient * np.asarray(gamma_i) * E_CHARGE * te)


def ion_momentum_flux(
    electron_density: float | np.ndarray,
    electron_temperature_ev: float | np.ndarray,
    ion_mass: float = DEUTERIUM_ION_MASS,
    ion_temperature_ev: float | np.ndarray = 0.0,
    ion_gamma: float = 1.0,
    ion_charge_state: float = 1.0,
    sheath_density_factor: float = 1.0,
) -> float | np.ndarray:
    '''Return the incident ion dynamic pressure m_i Gamma_i c_s.'''

    gamma_i = sheath_particle_flux(
        electron_density,
        electron_temperature_ev,
        ion_mass=ion_mass,
        ion_temperature_ev=ion_temperature_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=ion_charge_state,
        sheath_density_factor=sheath_density_factor,
    )
    sound_speed = bohm_sound_speed(
        electron_temperature_ev,
        ion_mass=ion_mass,
        ion_temperature_ev=ion_temperature_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=ion_charge_state,
    )
    return _maybe_scalar(ion_mass * np.asarray(gamma_i) * np.asarray(sound_speed))


def plasma_sheath_surface_loads(
    electron_density: float | np.ndarray,
    electron_temperature_ev: float | np.ndarray,
    ion_mass: float = DEUTERIUM_ION_MASS,
    ion_temperature_ev: float | np.ndarray = 0.0,
    ion_gamma: float = 1.0,
    ion_charge_state: float = 1.0,
    sheath_density_factor: float = 1.0,
    heat_transmission_coefficient: float = 7.0,
    net_current_density: float | np.ndarray = 0.0,
) -> PlasmaLiquidInterfaceLoads:
    '''Return reduced sheath loads to apply at a liquid-metal interface.'''

    particle_flux = np.asarray(
        sheath_particle_flux(
            electron_density,
            electron_temperature_ev,
            ion_mass=ion_mass,
            ion_temperature_ev=ion_temperature_ev,
            ion_gamma=ion_gamma,
            ion_charge_state=ion_charge_state,
            sheath_density_factor=sheath_density_factor,
        )
    )
    current_density = np.asarray(
        ion_saturation_current_density(
            electron_density,
            electron_temperature_ev,
            ion_mass=ion_mass,
            ion_temperature_ev=ion_temperature_ev,
            ion_gamma=ion_gamma,
            ion_charge_state=ion_charge_state,
            sheath_density_factor=sheath_density_factor,
        )
    )
    heat_flux = np.asarray(
        sheath_heat_transmission_flux(
            electron_density,
            electron_temperature_ev,
            ion_mass=ion_mass,
            ion_temperature_ev=ion_temperature_ev,
            ion_gamma=ion_gamma,
            ion_charge_state=ion_charge_state,
            sheath_density_factor=sheath_density_factor,
            heat_transmission_coefficient=heat_transmission_coefficient,
        )
    )
    normal_stress = np.asarray(
        ion_momentum_flux(
            electron_density,
            electron_temperature_ev,
            ion_mass=ion_mass,
            ion_temperature_ev=ion_temperature_ev,
            ion_gamma=ion_gamma,
            ion_charge_state=ion_charge_state,
            sheath_density_factor=sheath_density_factor,
        )
    )
    net_current = _finite_array("net_current_density", net_current_density)
    return PlasmaLiquidInterfaceLoads(
        particle_flux=particle_flux,
        incident_mass_flux=ion_mass * particle_flux,
        ion_current_density=current_density,
        net_current_density=np.broadcast_to(net_current, particle_flux.shape).copy(),
        heat_flux=heat_flux,
        normal_stress=normal_stress,
    )


def net_surface_heat_flux(
    plasma_heat_flux: float | np.ndarray,
    evaporation_mass_flux: float | np.ndarray = 0.0,
    latent_heat: float = 0.0,
) -> float | np.ndarray:
    '''Return heat flux retained by the liquid after evaporative cooling.

    `evaporation_mass_flux` is positive for material leaving the liquid surface.
    '''

    _require_nonnegative("latent_heat", latent_heat)
    heat = _finite_array("plasma_heat_flux", plasma_heat_flux)
    mdot = _nonnegative_array("evaporation_mass_flux", evaporation_mass_flux)
    return _maybe_scalar(heat - mdot * latent_heat)


def net_surface_mass_flux(
    incident_mass_flux: float | np.ndarray,
    evaporation_mass_flux: float | np.ndarray = 0.0,
) -> float | np.ndarray:
    '''Return positive-into-liquid mass flux after subtracting evaporation.'''

    incident = _nonnegative_array("incident_mass_flux", incident_mass_flux)
    evaporation = _nonnegative_array("evaporation_mass_flux", evaporation_mass_flux)
    return _maybe_scalar(incident - evaporation)


def surface_current_ohmic_heat_flux(
    current_density: float | np.ndarray,
    electrical_conductivity: float,
    conducting_thickness: float,
) -> float | np.ndarray:
    '''Return q_J = j^2 h / sigma for a uniform current across a layer.'''

    current = _finite_array("current_density", current_density)
    _require_positive("electrical_conductivity", electrical_conductivity)
    _require_positive("conducting_thickness", conducting_thickness)
    return _maybe_scalar(current**2 * conducting_thickness / electrical_conductivity)


def thermal_diffusivity(
    thermal_conductivity: float,
    density: float,
    specific_heat: float,
) -> float:
    '''Return thermal diffusivity alpha = k / (rho cp).'''

    _require_positive("thermal_conductivity", thermal_conductivity)
    _require_positive("density", density)
    _require_positive("specific_heat", specific_heat)
    return float(thermal_conductivity / (density * specific_heat))


def thermal_penetration_depth(
    time: float | np.ndarray,
    diffusivity: float,
    factor: float = 2.0,
) -> float | np.ndarray:
    '''Return a diffusion length factor*sqrt(alpha*t).'''

    t = _nonnegative_array("time", time)
    _require_positive("diffusivity", diffusivity)
    _require_positive("factor", factor)
    return _maybe_scalar(factor * np.sqrt(diffusivity * t))


def lumped_surface_temperature_rise(
    net_heat_flux: float | np.ndarray,
    time: float | np.ndarray,
    layer_thickness: float,
    density: float,
    specific_heat: float,
) -> float | np.ndarray:
    '''Return Delta T = q t / (rho cp h) for a well-mixed surface layer.'''

    heat = _finite_array("net_heat_flux", net_heat_flux)
    t = _nonnegative_array("time", time)
    _require_positive("layer_thickness", layer_thickness)
    _require_positive("density", density)
    _require_positive("specific_heat", specific_heat)
    return _maybe_scalar(heat * t / (density * specific_heat * layer_thickness))


def lumped_layer_energy_per_area(
    temperature_rise: float | np.ndarray,
    layer_thickness: float,
    density: float,
    specific_heat: float,
) -> float | np.ndarray:
    '''Return rho cp h Delta T, the sensible energy per unit surface area.'''

    delta_t = _finite_array("temperature_rise", temperature_rise)
    _require_positive("layer_thickness", layer_thickness)
    _require_positive("density", density)
    _require_positive("specific_heat", specific_heat)
    return _maybe_scalar(density * specific_heat * layer_thickness * delta_t)


def semi_infinite_constant_flux_temperature_rise(
    net_heat_flux: float | np.ndarray,
    time: float | np.ndarray,
    thermal_conductivity: float,
    density: float,
    specific_heat: float,
) -> float | np.ndarray:
    '''Return surface Delta T for a semi-infinite body with constant heat flux.

    For heat flux `q` applied at `x=0` into a semi-infinite liquid initially at
    uniform temperature, the analytic surface response is

        Delta T_s(t) = 2 q sqrt(alpha t / pi) / k,

    where `alpha = k/(rho cp)`. Positive `q` heats the liquid.
    '''

    heat = _finite_array("net_heat_flux", net_heat_flux)
    t = _nonnegative_array("time", time)
    _require_positive("thermal_conductivity", thermal_conductivity)
    alpha = thermal_diffusivity(thermal_conductivity, density, specific_heat)
    return _maybe_scalar(2.0 * heat * np.sqrt(alpha * t / np.pi) / thermal_conductivity)


def classify_thermal_response_regime(
    time: float | np.ndarray,
    layer_thickness: float,
    thermal_conductivity: float,
    density: float,
    specific_heat: float,
    semi_infinite_fraction: float = 0.25,
) -> str | np.ndarray:
    '''Classify thermal penetration relative to a liquid layer thickness.'''

    _require_positive("layer_thickness", layer_thickness)
    _require_positive("semi_infinite_fraction", semi_infinite_fraction)
    if semi_infinite_fraction >= 1.0:
        raise ValueError("semi_infinite_fraction must be less than 1")
    alpha = thermal_diffusivity(thermal_conductivity, density, specific_heat)
    depth = thermal_penetration_depth(time, alpha)
    labels = np.where(
        np.asarray(depth) < semi_infinite_fraction * layer_thickness,
        "semi-infinite",
        np.where(np.asarray(depth) < layer_thickness, "transient-slab", "through-thickness"),
    )
    if labels.ndim == 0:
        return str(labels.item())
    return labels


def marangoni_shear_stress(
    surface_temperature_gradient: float | np.ndarray,
    surface_tension_temperature_coefficient: float,
) -> float | np.ndarray:
    '''Return tau_s = (d sigma / dT) grad_s(T).'''

    gradient = _finite_array("surface_temperature_gradient", surface_temperature_gradient)
    coefficient = _finite_array(
        "surface_tension_temperature_coefficient",
        surface_tension_temperature_coefficient,
    )
    return _maybe_scalar(coefficient * gradient)


def thermocapillary_couette_profile(
    normal_coordinate: float | np.ndarray,
    layer_thickness: float,
    dynamic_viscosity: float,
    surface_temperature_gradient: float | np.ndarray,
    surface_tension_temperature_coefficient: float,
) -> float | np.ndarray:
    '''Return u(z) = tau_s z / mu for no-slip bottom and shear-driven surface.

    `normal_coordinate` is measured from the no-slip wall at `z=0` to the free
    surface at `z=layer_thickness`.
    '''

    z = _nonnegative_array("normal_coordinate", normal_coordinate)
    _require_positive("layer_thickness", layer_thickness)
    _require_positive("dynamic_viscosity", dynamic_viscosity)
    if np.any(z > layer_thickness):
        raise ValueError("normal_coordinate must not exceed layer_thickness")
    tau = marangoni_shear_stress(surface_temperature_gradient, surface_tension_temperature_coefficient)
    return _maybe_scalar(np.asarray(tau) * z / dynamic_viscosity)


def thermocapillary_couette_mean_velocity(
    layer_thickness: float,
    dynamic_viscosity: float,
    surface_temperature_gradient: float | np.ndarray,
    surface_tension_temperature_coefficient: float,
) -> float | np.ndarray:
    '''Return depth-mean velocity for the thermocapillary Couette profile.'''

    surface_velocity = thermocapillary_couette_profile(
        layer_thickness,
        layer_thickness,
        dynamic_viscosity,
        surface_temperature_gradient,
        surface_tension_temperature_coefficient,
    )
    return _maybe_scalar(0.5 * np.asarray(surface_velocity))


def capillary_gravity_stiffness(
    wavenumber: float | np.ndarray,
    density: float,
    surface_tension: float,
    gravity: float = STANDARD_GRAVITY,
    magnetic_stiffness: float | np.ndarray = 0.0,
) -> float | np.ndarray:
    '''Return static surface stiffness sigma k^2 + rho g + k_m.'''

    k = _positive_array("wavenumber", wavenumber)
    _require_positive("density", density)
    _require_positive("surface_tension", surface_tension)
    _require_nonnegative("gravity", gravity)
    magnetic = _nonnegative_array("magnetic_stiffness", magnetic_stiffness)
    return _maybe_scalar(surface_tension * k**2 + density * gravity + magnetic)


def static_capillary_gravity_displacement(
    normal_stress: float | np.ndarray,
    wavenumber: float | np.ndarray,
    density: float,
    surface_tension: float,
    gravity: float = STANDARD_GRAVITY,
    magnetic_stiffness: float | np.ndarray = 0.0,
) -> float | np.ndarray:
    '''Return eta = p_n / (sigma k^2 + rho g + k_m) for a linear surface mode.'''

    stress = _finite_array("normal_stress", normal_stress)
    stiffness = capillary_gravity_stiffness(
        wavenumber,
        density,
        surface_tension,
        gravity=gravity,
        magnetic_stiffness=magnetic_stiffness,
    )
    return _maybe_scalar(stress / np.asarray(stiffness))


def capillary_gravity_wave_frequency(
    wavenumber: float | np.ndarray,
    density: float,
    surface_tension: float,
    gravity: float = STANDARD_GRAVITY,
) -> float | np.ndarray:
    '''Return omega = sqrt(g k + sigma k^3/rho) for an inviscid deep surface mode.'''

    k = _positive_array("wavenumber", wavenumber)
    _require_positive("density", density)
    _require_positive("surface_tension", surface_tension)
    _require_nonnegative("gravity", gravity)
    return _maybe_scalar(np.sqrt(gravity * k + surface_tension * k**3 / density))


def magnetic_pressure_stiffness(
    wavenumber: float | np.ndarray,
    magnetic_field: float | np.ndarray,
    magnetic_permeability: float = VACUUM_PERMEABILITY,
    coefficient: float = 1.0,
) -> float | np.ndarray:
    '''Return an order-of-magnitude magnetic surface stiffness k B^2 / mu.

    This is a linear pressure-scale diagnostic for regime triage. The
    coefficient is intentionally explicit because the exact prefactor and sign
    depend on magnetic geometry and wall/free-surface boundary conditions.
    '''

    k = _positive_array("wavenumber", wavenumber)
    field = _finite_array("magnetic_field", magnetic_field)
    _require_positive("magnetic_permeability", magnetic_permeability)
    _require_nonnegative("coefficient", coefficient)
    return _maybe_scalar(coefficient * k * field**2 / magnetic_permeability)


def low_rm_magnetic_damping_rate(
    magnetic_field: float | np.ndarray,
    electrical_conductivity: float,
    density: float,
) -> float | np.ndarray:
    '''Return the low-Rm bulk Lorentz damping rate sigma B^2 / rho.'''

    field = _finite_array("magnetic_field", magnetic_field)
    _require_nonnegative("electrical_conductivity", electrical_conductivity)
    _require_positive("density", density)
    return _maybe_scalar(electrical_conductivity * field**2 / density)


def linear_surface_wave_frequency(
    wavenumber: float | np.ndarray,
    density: float,
    surface_tension: float,
    gravity: float = STANDARD_GRAVITY,
    magnetic_stiffness: float | np.ndarray = 0.0,
) -> float | np.ndarray:
    '''Return omega^2 = k (sigma k^2 + rho g + k_m) / rho.'''

    k = _positive_array("wavenumber", wavenumber)
    stiffness = capillary_gravity_stiffness(
        k,
        density,
        surface_tension,
        gravity=gravity,
        magnetic_stiffness=magnetic_stiffness,
    )
    return _maybe_scalar(np.sqrt(k * np.asarray(stiffness) / density))


def linear_surface_wave_period(
    wavenumber: float | np.ndarray,
    density: float,
    surface_tension: float,
    gravity: float = STANDARD_GRAVITY,
    magnetic_stiffness: float | np.ndarray = 0.0,
) -> float | np.ndarray:
    '''Return 2 pi / omega for the linear capillary-gravity-magnetic mode.'''

    omega = linear_surface_wave_frequency(
        wavenumber,
        density,
        surface_tension,
        gravity=gravity,
        magnetic_stiffness=magnetic_stiffness,
    )
    return _maybe_scalar(2.0 * np.pi / np.asarray(omega))


def surface_wave_explicit_timestep(
    angular_frequency: float | np.ndarray,
    damping_rate: float | np.ndarray = 0.0,
    safety_factor: float = 0.05,
) -> float | np.ndarray:
    '''Return a conservative explicit time-step gate safety/max(omega, gamma).'''

    omega = _positive_array("angular_frequency", angular_frequency)
    damping = _nonnegative_array("damping_rate", damping_rate)
    _require_positive("safety_factor", safety_factor)
    if safety_factor >= 1.0:
        raise ValueError("safety_factor must be less than 1")
    rate = np.maximum(omega, damping)
    return _maybe_scalar(safety_factor / rate)


def free_surface_linear_metrics(
    normal_stress: float | np.ndarray,
    wavenumber: float | np.ndarray,
    density: float,
    surface_tension: float,
    layer_thickness: float,
    exposure_time: float | np.ndarray,
    diffusivity: float,
    gravity: float = STANDARD_GRAVITY,
    magnetic_field: float | np.ndarray = 0.0,
    electrical_conductivity: float = 0.0,
    magnetic_permeability: float = VACUUM_PERMEABILITY,
    magnetic_stiffness_coefficient: float = 1.0,
    timestep_safety_factor: float = 0.05,
) -> FreeSurfaceLinearMetrics:
    '''Return linear response, stiffness, damping, and time-step gate metrics.

    `normal_stress` is positive into the liquid and produces a signed static
    displacement under the caller's surface-normal convention. The magnetic
    stiffness term is an order-of-magnitude pressure scale; set
    `magnetic_stiffness_coefficient=0` to disable it.
    '''

    stress = _finite_array("normal_stress", normal_stress)
    k = _positive_array("wavenumber", wavenumber)
    _require_positive("layer_thickness", layer_thickness)
    t = _nonnegative_array("exposure_time", exposure_time)
    _require_positive("diffusivity", diffusivity)
    magnetic_stiffness = magnetic_pressure_stiffness(
        k,
        magnetic_field,
        magnetic_permeability=magnetic_permeability,
        coefficient=magnetic_stiffness_coefficient,
    )
    stiffness = capillary_gravity_stiffness(
        k,
        density,
        surface_tension,
        gravity=gravity,
        magnetic_stiffness=magnetic_stiffness,
    )
    omega = linear_surface_wave_frequency(
        k,
        density,
        surface_tension,
        gravity=gravity,
        magnetic_stiffness=magnetic_stiffness,
    )
    period = 2.0 * np.pi / np.asarray(omega)
    damping = low_rm_magnetic_damping_rate(magnetic_field, electrical_conductivity, density)
    damping_ratio = np.asarray(damping) / (2.0 * np.asarray(omega))
    explicit_dt = surface_wave_explicit_timestep(
        omega,
        damping_rate=damping,
        safety_factor=timestep_safety_factor,
    )
    displacement = stress / np.asarray(stiffness)
    thermal_fraction = thermal_penetration_depth(t, diffusivity) / layer_thickness
    return FreeSurfaceLinearMetrics(
        wavenumber=np.asarray(k),
        wavelength=2.0 * np.pi / np.asarray(k),
        stiffness=np.asarray(stiffness),
        magnetic_stiffness=np.asarray(magnetic_stiffness),
        angular_frequency=np.asarray(omega),
        period=np.asarray(period),
        magnetic_damping_rate=np.asarray(damping),
        magnetic_damping_ratio=np.asarray(damping_ratio),
        explicit_timestep=np.asarray(explicit_dt),
        static_displacement=np.asarray(displacement),
        displacement_fraction=np.abs(np.asarray(displacement)) / layer_thickness,
        thermal_penetration_fraction=np.asarray(thermal_fraction),
    )


def _finite_array(name: str, values: float | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _positive_array(name: str, values: float | np.ndarray) -> np.ndarray:
    arr = _finite_array(name, values)
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive")
    return arr


def _nonnegative_array(name: str, values: float | np.ndarray) -> np.ndarray:
    arr = _finite_array(name, values)
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative")
    return arr


def _require_positive(name: str, value: float) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive")


def _require_nonnegative(name: str, value: float) -> None:
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be non-negative")


def _maybe_scalar(values: np.ndarray) -> float | np.ndarray:
    if values.ndim == 0:
        return float(values)
    return values
