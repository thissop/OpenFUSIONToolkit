import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        DEUTERIUM_ION_MASS,
        E_CHARGE,
        STANDARD_GRAVITY,
        VACUUM_PERMEABILITY,
        bohm_sound_speed,
        capillary_gravity_stiffness,
        capillary_gravity_wave_frequency,
        classify_thermal_response_regime,
        free_surface_linear_metrics,
        ion_momentum_flux,
        ion_saturation_current_density,
        linear_surface_wave_frequency,
        linear_surface_wave_period,
        lumped_layer_energy_per_area,
        lumped_surface_temperature_rise,
        low_rm_magnetic_damping_rate,
        magnetic_pressure_stiffness,
        marangoni_shear_stress,
        net_surface_heat_flux,
        net_surface_mass_flux,
        plasma_sheath_surface_loads,
        representative_liquid_metal_thermal_properties,
        semi_infinite_constant_flux_temperature_rise,
        sheath_heat_transmission_flux,
        sheath_particle_flux,
        static_capillary_gravity_displacement,
        surface_current_ohmic_heat_flux,
        surface_wave_explicit_timestep,
        thermocapillary_couette_mean_velocity,
        thermocapillary_couette_profile,
        thermal_diffusivity,
        thermal_penetration_depth,
    )
except FileNotFoundError:
    # Source-tree fallback for Mac-side development shells where the compiled
    # OFT shared library is not present next to the Python package.
    import importlib.util
    import sys
    from pathlib import Path

    interface_path = (
        Path(__file__).resolve().parents[2] / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "plasma_interface.py"
    )
    spec = importlib.util.spec_from_file_location("lm_mhd_plasma_interface", interface_path)
    plasma_interface = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = plasma_interface
    spec.loader.exec_module(plasma_interface)
    DEUTERIUM_ION_MASS = plasma_interface.DEUTERIUM_ION_MASS
    E_CHARGE = plasma_interface.E_CHARGE
    STANDARD_GRAVITY = plasma_interface.STANDARD_GRAVITY
    VACUUM_PERMEABILITY = plasma_interface.VACUUM_PERMEABILITY
    bohm_sound_speed = plasma_interface.bohm_sound_speed
    capillary_gravity_stiffness = plasma_interface.capillary_gravity_stiffness
    capillary_gravity_wave_frequency = plasma_interface.capillary_gravity_wave_frequency
    classify_thermal_response_regime = plasma_interface.classify_thermal_response_regime
    free_surface_linear_metrics = plasma_interface.free_surface_linear_metrics
    ion_momentum_flux = plasma_interface.ion_momentum_flux
    ion_saturation_current_density = plasma_interface.ion_saturation_current_density
    linear_surface_wave_frequency = plasma_interface.linear_surface_wave_frequency
    linear_surface_wave_period = plasma_interface.linear_surface_wave_period
    lumped_layer_energy_per_area = plasma_interface.lumped_layer_energy_per_area
    lumped_surface_temperature_rise = plasma_interface.lumped_surface_temperature_rise
    low_rm_magnetic_damping_rate = plasma_interface.low_rm_magnetic_damping_rate
    magnetic_pressure_stiffness = plasma_interface.magnetic_pressure_stiffness
    marangoni_shear_stress = plasma_interface.marangoni_shear_stress
    net_surface_heat_flux = plasma_interface.net_surface_heat_flux
    net_surface_mass_flux = plasma_interface.net_surface_mass_flux
    plasma_sheath_surface_loads = plasma_interface.plasma_sheath_surface_loads
    representative_liquid_metal_thermal_properties = plasma_interface.representative_liquid_metal_thermal_properties
    semi_infinite_constant_flux_temperature_rise = plasma_interface.semi_infinite_constant_flux_temperature_rise
    sheath_heat_transmission_flux = plasma_interface.sheath_heat_transmission_flux
    sheath_particle_flux = plasma_interface.sheath_particle_flux
    static_capillary_gravity_displacement = plasma_interface.static_capillary_gravity_displacement
    surface_current_ohmic_heat_flux = plasma_interface.surface_current_ohmic_heat_flux
    surface_wave_explicit_timestep = plasma_interface.surface_wave_explicit_timestep
    thermocapillary_couette_mean_velocity = plasma_interface.thermocapillary_couette_mean_velocity
    thermocapillary_couette_profile = plasma_interface.thermocapillary_couette_profile
    thermal_diffusivity = plasma_interface.thermal_diffusivity
    thermal_penetration_depth = plasma_interface.thermal_penetration_depth


def test_bohm_particle_current_and_momentum_identities():
    density = 2.0e19
    te_ev = 25.0
    ti_ev = 5.0
    ion_gamma = 3.0
    z_ion = 1.0

    sound_speed = bohm_sound_speed(
        te_ev,
        ion_temperature_ev=ti_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=z_ion,
    )
    expected_cs = np.sqrt(E_CHARGE * (z_ion * te_ev + ion_gamma * ti_ev) / DEUTERIUM_ION_MASS)
    np.testing.assert_allclose(sound_speed, expected_cs)

    particle_flux = sheath_particle_flux(
        density,
        te_ev,
        ion_temperature_ev=ti_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=z_ion,
    )
    current = ion_saturation_current_density(
        density,
        te_ev,
        ion_temperature_ev=ti_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=z_ion,
    )
    momentum = ion_momentum_flux(
        density,
        te_ev,
        ion_temperature_ev=ti_ev,
        ion_gamma=ion_gamma,
        ion_charge_state=z_ion,
    )

    np.testing.assert_allclose(particle_flux, density * sound_speed)
    np.testing.assert_allclose(current, z_ion * E_CHARGE * particle_flux)
    np.testing.assert_allclose(momentum, DEUTERIUM_ION_MASS * particle_flux * sound_speed)


def test_sheath_surface_loads_default_to_floating_net_current():
    density = np.array([1.0e19, 2.0e19])
    te_ev = np.array([10.0, 30.0])
    gamma_sh = 6.5
    imposed_current = np.array([0.0, 25.0])

    floating = plasma_sheath_surface_loads(density, te_ev, heat_transmission_coefficient=gamma_sh)
    biased = plasma_sheath_surface_loads(
        density,
        te_ev,
        heat_transmission_coefficient=gamma_sh,
        net_current_density=imposed_current,
    )

    np.testing.assert_allclose(floating.net_current_density, 0.0)
    np.testing.assert_allclose(biased.net_current_density, imposed_current)
    np.testing.assert_allclose(floating.ion_current_density, E_CHARGE * floating.particle_flux)
    np.testing.assert_allclose(floating.incident_mass_flux, DEUTERIUM_ION_MASS * floating.particle_flux)
    np.testing.assert_allclose(floating.heat_flux, gamma_sh * floating.particle_flux * E_CHARGE * te_ev)
    assert np.all(floating.normal_stress > 0.0)


def test_sheath_heat_flux_is_transmission_coefficient_times_particle_energy():
    density = 1.5e19
    te_ev = 18.0
    gamma_sh = 7.2

    q = sheath_heat_transmission_flux(density, te_ev, heat_transmission_coefficient=gamma_sh)
    gamma_i = sheath_particle_flux(density, te_ev)

    assert q == pytest.approx(gamma_sh * gamma_i * E_CHARGE * te_ev)


def test_evaporation_energy_and_mass_accounting_use_positive_outflow():
    plasma_heat = np.array([2.0e6, 5.0e5])
    evaporation = np.array([1.0e-4, 4.0e-4])
    latent_heat = 2.0e7

    retained_heat = net_surface_heat_flux(plasma_heat, evaporation, latent_heat)
    net_mass = net_surface_mass_flux(np.array([3.0e-4, 3.0e-4]), evaporation)

    np.testing.assert_allclose(retained_heat, plasma_heat - evaporation * latent_heat)
    np.testing.assert_allclose(net_mass, np.array([2.0e-4, -1.0e-4]))
    with pytest.raises(ValueError, match="evaporation_mass_flux"):
        net_surface_heat_flux(plasma_heat, -1.0e-4, latent_heat)


def test_surface_current_ohmic_heat_flux_is_positive_and_even_in_current():
    current = np.array([-2.0e4, 0.0, 3.0e4])
    sigma = 8.0e5
    thickness = 0.02

    heat = surface_current_ohmic_heat_flux(current, sigma, thickness)

    np.testing.assert_allclose(heat, current**2 * thickness / sigma)
    assert np.all(heat >= 0.0)
    with pytest.raises(ValueError, match="electrical_conductivity"):
        surface_current_ohmic_heat_flux(current, 0.0, thickness)


def test_static_capillary_gravity_response_and_wave_frequency():
    k = np.array([5.0, 10.0])
    density = 500.0
    surface_tension = 0.38
    magnetic_stiffness = np.array([100.0, 200.0])
    stress = np.array([20.0, 40.0])

    stiffness = capillary_gravity_stiffness(
        k,
        density,
        surface_tension,
        magnetic_stiffness=magnetic_stiffness,
    )
    displacement = static_capillary_gravity_displacement(
        stress,
        k,
        density,
        surface_tension,
        magnetic_stiffness=magnetic_stiffness,
    )
    omega = capillary_gravity_wave_frequency(k, density, surface_tension)

    np.testing.assert_allclose(stiffness, surface_tension * k**2 + density * STANDARD_GRAVITY + magnetic_stiffness)
    np.testing.assert_allclose(displacement, stress / stiffness)
    np.testing.assert_allclose(omega**2, STANDARD_GRAVITY * k + surface_tension * k**3 / density)


def test_linear_surface_wave_and_magnetic_scales_are_consistent():
    k = np.array([20.0, 80.0])
    density = 500.0
    surface_tension = 0.38
    magnetic_field = 2.0
    conductivity = 3.2e6
    safety = 0.05

    magnetic_stiffness = magnetic_pressure_stiffness(k, magnetic_field)
    stiffness = capillary_gravity_stiffness(
        k,
        density,
        surface_tension,
        magnetic_stiffness=magnetic_stiffness,
    )
    omega = linear_surface_wave_frequency(
        k,
        density,
        surface_tension,
        magnetic_stiffness=magnetic_stiffness,
    )
    period = linear_surface_wave_period(
        k,
        density,
        surface_tension,
        magnetic_stiffness=magnetic_stiffness,
    )
    damping = low_rm_magnetic_damping_rate(magnetic_field, conductivity, density)
    timestep = surface_wave_explicit_timestep(omega, damping_rate=damping, safety_factor=safety)

    np.testing.assert_allclose(magnetic_stiffness, k * magnetic_field**2 / VACUUM_PERMEABILITY)
    np.testing.assert_allclose(omega**2, k * stiffness / density)
    np.testing.assert_allclose(period, 2.0 * np.pi / omega)
    np.testing.assert_allclose(damping, conductivity * magnetic_field**2 / density)
    np.testing.assert_allclose(timestep, safety / np.maximum(omega, damping))


def test_free_surface_linear_metrics_close_response_identities():
    props = representative_liquid_metal_thermal_properties()["Li"]
    k = np.array([30.0, 60.0])
    stress = np.array([12.0, 18.0])
    layer = 0.003
    exposure = np.array([0.002, 0.008])
    magnetic_field = 1.5
    conductivity = 3.2e6

    metrics = free_surface_linear_metrics(
        stress,
        k,
        props.density,
        props.surface_tension,
        layer,
        exposure,
        props.thermal_diffusivity,
        magnetic_field=magnetic_field,
        electrical_conductivity=conductivity,
    )

    np.testing.assert_allclose(metrics.wavelength, 2.0 * np.pi / k)
    np.testing.assert_allclose(metrics.angular_frequency**2, k * metrics.stiffness / props.density)
    np.testing.assert_allclose(metrics.static_displacement, stress / metrics.stiffness)
    np.testing.assert_allclose(metrics.displacement_fraction, np.abs(metrics.static_displacement) / layer)
    np.testing.assert_allclose(
        metrics.thermal_penetration_fraction,
        thermal_penetration_depth(exposure, props.thermal_diffusivity) / layer,
    )
    np.testing.assert_allclose(
        metrics.magnetic_damping_ratio,
        metrics.magnetic_damping_rate / (2.0 * metrics.angular_frequency),
    )
    with pytest.raises(ValueError, match="safety_factor"):
        surface_wave_explicit_timestep(metrics.angular_frequency, safety_factor=1.0)


def test_interface_helpers_reject_nonphysical_inputs():
    with pytest.raises(ValueError, match="electron_temperature_ev"):
        bohm_sound_speed(0.0)
    with pytest.raises(ValueError, match="electron_density"):
        sheath_particle_flux(-1.0, 10.0)
    with pytest.raises(ValueError, match="wavenumber"):
        capillary_gravity_stiffness(0.0, density=500.0, surface_tension=0.38)
    with pytest.raises(ValueError, match="latent_heat"):
        net_surface_heat_flux(1.0e6, evaporation_mass_flux=1.0e-4, latent_heat=-1.0)


def test_lumped_surface_temperature_rise_closes_energy_per_area():
    q = np.array([2.0e6, -5.0e5])
    time = np.array([0.02, 0.04])
    thickness = 2.5e-3
    density = 500.0
    cp = 4200.0

    delta_t = lumped_surface_temperature_rise(q, time, thickness, density, cp)
    energy = lumped_layer_energy_per_area(delta_t, thickness, density, cp)

    np.testing.assert_allclose(energy, q * time)


def test_semi_infinite_constant_flux_response_has_sqrt_time_scaling():
    q = 1.2e6
    density = 500.0
    cp = 4200.0
    conductivity = 50.0
    time = np.array([0.0, 0.01, 0.04])
    alpha = thermal_diffusivity(conductivity, density, cp)

    delta_t = semi_infinite_constant_flux_temperature_rise(q, time, conductivity, density, cp)

    np.testing.assert_allclose(alpha, conductivity / (density * cp))
    np.testing.assert_allclose(delta_t[0], 0.0)
    np.testing.assert_allclose(delta_t[2], 2.0 * delta_t[1])
    expected = 2.0 * q * np.sqrt(alpha * time / np.pi) / conductivity
    np.testing.assert_allclose(delta_t, expected)


def test_thermal_penetration_depth_and_regime_classifier():
    props = representative_liquid_metal_thermal_properties()["Li"]
    alpha = props.thermal_diffusivity
    layer = 0.01
    times = np.array(
        [
            (0.1 * layer / 2.0) ** 2 / alpha,
            (0.5 * layer / 2.0) ** 2 / alpha,
            (1.5 * layer / 2.0) ** 2 / alpha,
        ]
    )

    depth = thermal_penetration_depth(times, alpha)
    labels = classify_thermal_response_regime(
        times,
        layer,
        props.thermal_conductivity,
        props.density,
        props.specific_heat,
    )

    np.testing.assert_allclose(depth / layer, np.array([0.1, 0.5, 1.5]))
    assert labels.tolist() == ["semi-infinite", "transient-slab", "through-thickness"]
    with pytest.raises(ValueError, match="semi_infinite_fraction"):
        classify_thermal_response_regime(times, layer, props.thermal_conductivity, props.density, props.specific_heat, 1.0)


def test_thermocapillary_couette_profile_is_linear_with_correct_sign():
    z = np.linspace(0.0, 0.004, 5)
    layer = 0.004
    viscosity = 5.0e-4
    d_sigma_d_t = -1.0e-4
    temperature_gradient = 2.5e3
    tau = marangoni_shear_stress(temperature_gradient, d_sigma_d_t)

    velocity = thermocapillary_couette_profile(z, layer, viscosity, temperature_gradient, d_sigma_d_t)

    np.testing.assert_allclose(tau, d_sigma_d_t * temperature_gradient)
    np.testing.assert_allclose(velocity, tau * z / viscosity)
    np.testing.assert_allclose(velocity[0], 0.0)
    assert velocity[-1] < 0.0
    with pytest.raises(ValueError, match="normal_coordinate"):
        thermocapillary_couette_profile(layer * 1.1, layer, viscosity, temperature_gradient, d_sigma_d_t)


def test_thermocapillary_couette_mean_is_half_surface_velocity():
    layer = 0.003
    viscosity = 4.5e-4
    d_sigma_d_t = -8.0e-5
    temperature_gradient = np.array([1.0e3, 2.0e3])

    surface_velocity = thermocapillary_couette_profile(
        layer,
        layer,
        viscosity,
        temperature_gradient,
        d_sigma_d_t,
    )
    mean_velocity = thermocapillary_couette_mean_velocity(layer, viscosity, temperature_gradient, d_sigma_d_t)

    np.testing.assert_allclose(mean_velocity, 0.5 * surface_velocity)
