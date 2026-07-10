import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        axisymmetric_boundary_current_integrals,
        axisymmetric_cell_volumes,
        axisymmetric_charge_conservation_residual,
        axisymmetric_current_from_flux_function,
        axisymmetric_electric_field_from_potential,
        axisymmetric_integrated_power,
        axisymmetric_integrated_torque,
        axisymmetric_net_boundary_current,
        axisymmetric_torque_density,
        axisymmetric_volume_integral,
        axisymmetric_wall_normal_current_extrema,
        charge_conservation_residual_2d,
        divergence_free_current_from_streamfunction,
        electric_power_density,
        insulating_wall_normal_gradient,
        insulating_wall_normal_flux,
        inductionless_power_balance_residual,
        joule_heating_density,
        lorentz_force_density,
        mechanical_power_density,
        motional_electric_field,
        ohms_law_current_density,
        potential_source_from_motional_emf,
        reconstruct_axisymmetric_inductionless_current,
        reconstruct_inductionless_current_2d,
        solve_potential_dirichlet_2d,
        solve_potential_neumann_2d,
        solve_variable_conductivity_dirichlet_2d,
        solve_variable_conductivity_neumann_2d,
        structured_gradient_2d,
        toroidal_electric_field_from_loop_voltage,
        wall_normal_current_extrema,
        weighted_potential_source_from_motional_emf,
    )
except FileNotFoundError:
    # Source-tree fallback for Mac-side development shells where the compiled
    # OFT shared library is not present next to the Python package.
    import importlib.util
    import sys
    from pathlib import Path

    inductionless_path = (
        Path(__file__).resolve().parents[2] / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "inductionless.py"
    )
    spec = importlib.util.spec_from_file_location("lm_mhd_inductionless", inductionless_path)
    inductionless = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = inductionless
    spec.loader.exec_module(inductionless)
    axisymmetric_cell_volumes = inductionless.axisymmetric_cell_volumes
    axisymmetric_boundary_current_integrals = inductionless.axisymmetric_boundary_current_integrals
    axisymmetric_charge_conservation_residual = inductionless.axisymmetric_charge_conservation_residual
    axisymmetric_current_from_flux_function = inductionless.axisymmetric_current_from_flux_function
    axisymmetric_electric_field_from_potential = inductionless.axisymmetric_electric_field_from_potential
    axisymmetric_integrated_power = inductionless.axisymmetric_integrated_power
    axisymmetric_integrated_torque = inductionless.axisymmetric_integrated_torque
    axisymmetric_net_boundary_current = inductionless.axisymmetric_net_boundary_current
    axisymmetric_torque_density = inductionless.axisymmetric_torque_density
    axisymmetric_volume_integral = inductionless.axisymmetric_volume_integral
    axisymmetric_wall_normal_current_extrema = inductionless.axisymmetric_wall_normal_current_extrema
    charge_conservation_residual_2d = inductionless.charge_conservation_residual_2d
    divergence_free_current_from_streamfunction = inductionless.divergence_free_current_from_streamfunction
    electric_power_density = inductionless.electric_power_density
    insulating_wall_normal_gradient = inductionless.insulating_wall_normal_gradient
    insulating_wall_normal_flux = inductionless.insulating_wall_normal_flux
    inductionless_power_balance_residual = inductionless.inductionless_power_balance_residual
    joule_heating_density = inductionless.joule_heating_density
    lorentz_force_density = inductionless.lorentz_force_density
    mechanical_power_density = inductionless.mechanical_power_density
    motional_electric_field = inductionless.motional_electric_field
    ohms_law_current_density = inductionless.ohms_law_current_density
    potential_source_from_motional_emf = inductionless.potential_source_from_motional_emf
    reconstruct_axisymmetric_inductionless_current = inductionless.reconstruct_axisymmetric_inductionless_current
    reconstruct_inductionless_current_2d = inductionless.reconstruct_inductionless_current_2d
    solve_potential_dirichlet_2d = inductionless.solve_potential_dirichlet_2d
    solve_potential_neumann_2d = inductionless.solve_potential_neumann_2d
    solve_variable_conductivity_dirichlet_2d = inductionless.solve_variable_conductivity_dirichlet_2d
    solve_variable_conductivity_neumann_2d = inductionless.solve_variable_conductivity_neumann_2d
    structured_gradient_2d = inductionless.structured_gradient_2d
    toroidal_electric_field_from_loop_voltage = inductionless.toroidal_electric_field_from_loop_voltage
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema
    weighted_potential_source_from_motional_emf = inductionless.weighted_potential_source_from_motional_emf


def test_ohms_law_cancels_motional_emf_with_matching_potential_gradient():
    velocity = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    magnetic_field = np.array([0.0, 0.0, 3.0])
    grad_phi = np.cross(velocity, magnetic_field)

    current = ohms_law_current_density(2.5, velocity, magnetic_field, grad_phi)

    np.testing.assert_allclose(current, 0.0, atol=1.0e-14)


def test_ohms_law_broadcasts_scalar_and_spatial_sigma():
    velocity = np.array([[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]])
    magnetic_field = np.array([0.0, 1.0, 0.0])
    sigma = np.array([[2.0, 3.0]])

    current = ohms_law_current_density(sigma, velocity, magnetic_field)

    expected = sigma[..., None] * np.cross(velocity, magnetic_field)
    np.testing.assert_allclose(current, expected)


def test_lorentz_force_is_cross_product_and_perpendicular():
    current = np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, 2.0]])
    magnetic_field = np.array([0.2, -0.3, 0.7])

    force = lorentz_force_density(current, magnetic_field)

    np.testing.assert_allclose(force, np.cross(current, magnetic_field))
    np.testing.assert_allclose(np.sum(force * current, axis=-1), 0.0, atol=1.0e-14)
    np.testing.assert_allclose(np.sum(force * magnetic_field, axis=-1), 0.0, atol=1.0e-14)


def test_joule_heating_is_nonnegative_and_rejects_bad_sigma():
    current = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 2.0]])

    heating = joule_heating_density(current, sigma=2.0)

    np.testing.assert_allclose(heating, np.array([12.5, 2.0]))
    assert np.all(heating >= 0.0)
    with pytest.raises(ValueError, match="sigma"):
        joule_heating_density(current, sigma=0.0)


def test_streamfunction_current_is_discretely_charge_conserving():
    x = np.linspace(0.0, 1.0, 65)
    z = np.linspace(0.0, 1.0, 65)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    streamfunction = np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)

    current = divergence_free_current_from_streamfunction(streamfunction, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)

    assert np.max(np.abs(residual[2:-2, 2:-2])) < 1.0e-11
    assert max(wall_current.values()) < 1.0e-12


def test_axisymmetric_flux_function_current_is_charge_conserving():
    r = np.linspace(1.0, 2.0, 65)
    z = np.linspace(-0.4, 0.4, 67)
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    flux = (r_grid - 1.0) ** 2 * (2.0 - r_grid) ** 2 * (z_grid + 0.4) ** 2 * (0.4 - z_grid) ** 2

    current = axisymmetric_current_from_flux_function(flux, r, z)
    residual = axisymmetric_charge_conservation_residual(current, r, z)
    wall_current = axisymmetric_wall_normal_current_extrema(current)
    boundary_current = axisymmetric_boundary_current_integrals(current, r, z)

    assert np.max(np.abs(residual[2:-2, 2:-2])) < 1.0e-10
    assert max(wall_current.values()) < 1.0e-12
    assert max(abs(value) for value in boundary_current.values()) < 1.0e-12
    assert abs(axisymmetric_net_boundary_current(current, r, z)) < 1.0e-12


def test_axisymmetric_boundary_current_integrals_use_outward_surface_metrics():
    r = np.linspace(1.0, 2.0, 17)
    z = np.linspace(-0.5, 0.5, 19)
    r_grid, _ = np.meshgrid(r, z, indexing="xy")
    radial_strength = 7.0
    vertical_strength = -2.0
    current = np.zeros(r_grid.shape + (3,))
    current[..., 0] = radial_strength / r_grid
    current[..., 2] = vertical_strength

    residual = axisymmetric_charge_conservation_residual(current, r, z)
    boundary_current = axisymmetric_boundary_current_integrals(current, r, z)

    height = z[-1] - z[0]
    annulus_area = np.pi * (r[-1] ** 2 - r[0] ** 2)
    expected_radial = 2.0 * np.pi * radial_strength * height
    expected_vertical = vertical_strength * annulus_area
    np.testing.assert_allclose(boundary_current["r_min"], -expected_radial)
    np.testing.assert_allclose(boundary_current["r_max"], expected_radial)
    np.testing.assert_allclose(boundary_current["z_min"], -expected_vertical)
    np.testing.assert_allclose(boundary_current["z_max"], expected_vertical)
    np.testing.assert_allclose(axisymmetric_net_boundary_current(current, r, z), 0.0, atol=1.0e-13)
    assert np.max(np.abs(residual)) < 1.0e-12
    with pytest.raises(ValueError, match="shape"):
        axisymmetric_boundary_current_integrals(current[:, :-1, :], r, z)


def test_axisymmetric_toroidal_electric_field_drives_jphi_without_divergence():
    r = np.linspace(1.0, 2.0, 17)
    z = np.linspace(-0.5, 0.5, 19)
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    sigma = 1.0 + 0.2 * (r_grid - 1.0)
    potential = np.zeros_like(r_grid)
    velocity = np.zeros(r_grid.shape + (3,))
    magnetic_field = np.zeros_like(velocity)
    magnetic_field[..., 2] = 0.3
    e_phi = 0.04

    current = reconstruct_axisymmetric_inductionless_current(
        potential,
        velocity,
        magnetic_field,
        sigma,
        r,
        z,
        toroidal_electric_field=e_phi,
    )
    residual = axisymmetric_charge_conservation_residual(current, r, z)
    force = lorentz_force_density(current, magnetic_field)

    np.testing.assert_allclose(current[..., 1], sigma * e_phi)
    np.testing.assert_allclose(current[..., 0], 0.0, atol=1.0e-14)
    np.testing.assert_allclose(current[..., 2], 0.0, atol=1.0e-14)
    np.testing.assert_allclose(residual, 0.0, atol=1.0e-14)
    np.testing.assert_allclose(force[..., 0], sigma * e_phi * magnetic_field[..., 2])


def test_axisymmetric_electric_field_from_potential_uses_minus_gradient_and_ephi():
    r = np.linspace(1.0, 2.0, 17)
    z = np.linspace(-0.5, 0.5, 19)
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    potential = 0.25 * r_grid**2 - 0.3 * z_grid
    e_phi = 0.04 * r_grid

    electric_field = axisymmetric_electric_field_from_potential(potential, r, z, toroidal_electric_field=e_phi)

    np.testing.assert_allclose(electric_field[..., 0], -0.5 * r_grid, atol=1.0e-14)
    np.testing.assert_allclose(electric_field[..., 1], e_phi)
    np.testing.assert_allclose(electric_field[..., 2], 0.3, atol=1.0e-14)


def test_inductionless_power_balance_identity_holds_pointwise():
    velocity = np.array(
        [
            [[0.2, -0.1, 0.4], [0.1, 0.3, -0.2]],
            [[-0.3, 0.2, 0.1], [0.4, -0.2, 0.05]],
        ]
    )
    magnetic_field = np.array([0.7, -0.4, 1.2])
    electric_field = np.array(
        [
            [[0.01, 0.03, -0.02], [0.02, -0.01, 0.04]],
            [[-0.03, 0.02, 0.01], [0.05, -0.02, -0.01]],
        ]
    )
    sigma = np.array([[2.0, 3.0], [4.0, 5.0]])
    current = sigma[..., None] * (electric_field + np.cross(velocity, magnetic_field))
    force = lorentz_force_density(current, magnetic_field)

    residual = inductionless_power_balance_residual(current, velocity, magnetic_field, sigma, electric_field)
    joule = joule_heating_density(current, sigma)
    mechanical = mechanical_power_density(force, velocity)
    electric = electric_power_density(current, electric_field)

    np.testing.assert_allclose(residual, 0.0, atol=1.0e-14)
    np.testing.assert_allclose(joule + mechanical, electric, atol=1.0e-14)
    with pytest.raises(ValueError, match="broadcast"):
        inductionless_power_balance_residual(current, np.zeros((3, 2, 3)), magnetic_field, sigma, electric_field)


def test_loop_voltage_sets_one_over_r_toroidal_electric_field():
    r = np.array([1.0, 1.5, 2.0])
    loop_voltage = 0.12

    e_phi = toroidal_electric_field_from_loop_voltage(loop_voltage, r)

    np.testing.assert_allclose(2.0 * np.pi * r * e_phi, loop_voltage)
    with pytest.raises(ValueError, match="r"):
        toroidal_electric_field_from_loop_voltage(loop_voltage, np.array([0.0, 1.0]))


def test_axisymmetric_cell_volumes_integrate_annulus_area():
    r = np.linspace(1.0, 2.0, 65)
    z = np.linspace(-0.5, 0.5, 33)

    volumes = axisymmetric_cell_volumes(r, z)

    expected = np.pi * (r[-1] ** 2 - r[0] ** 2) * (z[-1] - z[0])
    np.testing.assert_allclose(np.sum(volumes), expected, rtol=1.0e-12, atol=0.0)


def test_mechanical_power_density_is_vector_dot_product():
    force = np.array(
        [
            [[1.0, 2.0, -1.0], [0.5, -2.0, 3.0]],
            [[-4.0, 0.0, 2.0], [1.5, 1.0, 0.0]],
        ]
    )
    velocity = np.array([2.0, -1.0, 0.5])

    power = mechanical_power_density(force, velocity)

    np.testing.assert_allclose(power, np.sum(force * velocity, axis=-1))
    with pytest.raises(ValueError, match="broadcast"):
        mechanical_power_density(force, np.zeros((3, 2, 3)))


def test_axisymmetric_volume_integral_and_integrated_power_use_annular_weights():
    r = np.linspace(1.0, 2.0, 65)
    z = np.linspace(-0.5, 0.5, 33)
    density = np.full((z.size, r.size), 2.5)

    expected_volume = np.pi * (r[-1] ** 2 - r[0] ** 2) * (z[-1] - z[0])

    np.testing.assert_allclose(axisymmetric_volume_integral(density, r, z), 2.5 * expected_volume)
    np.testing.assert_allclose(axisymmetric_integrated_power(density, r, z), 2.5 * expected_volume)
    with pytest.raises(ValueError, match="shape"):
        axisymmetric_volume_integral(density[:, :-1], r, z)


def test_axisymmetric_torque_density_and_integral_are_about_symmetry_axis():
    r = np.linspace(1.0, 2.0, 65)
    z = np.linspace(-0.5, 0.5, 33)
    r_grid, _ = np.meshgrid(r, z, indexing="xy")
    force = np.zeros(r_grid.shape + (3,))
    force[..., 1] = 3.0 / r_grid

    torque_density = axisymmetric_torque_density(force, r, z)
    torque = axisymmetric_integrated_torque(force, r, z)

    expected_volume = np.pi * (r[-1] ** 2 - r[0] ** 2) * (z[-1] - z[0])
    np.testing.assert_allclose(torque_density, 3.0)
    np.testing.assert_allclose(torque, 3.0 * expected_volume)
    with pytest.raises(ValueError, match="shape"):
        axisymmetric_torque_density(force[:, :-1, :], r, z)


def test_structured_residual_rejects_shape_mismatch():
    current = np.zeros((4, 5, 3))
    with pytest.raises(ValueError, match="shape"):
        charge_conservation_residual_2d(current, np.linspace(0.0, 1.0, 6), np.linspace(0.0, 1.0, 4))


def test_dirichlet_poisson_solver_recovers_sine_mms():
    x = np.linspace(0.0, 1.0, 51)
    z = np.linspace(0.0, 1.0, 49)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    source = -2.0 * np.pi**2 * phi_exact

    phi = solve_potential_dirichlet_2d(source, x, z, boundary_values=0.0)

    max_error = np.max(np.abs(phi - phi_exact))
    rms_error = np.sqrt(np.mean((phi - phi_exact) ** 2))
    assert max_error < 4.0e-4
    assert rms_error < 2.0e-4


def test_dirichlet_poisson_solver_rejects_nonuniform_grid():
    source = np.zeros((5, 5))
    x = np.array([0.0, 0.2, 0.5, 0.75, 1.0])
    z = np.linspace(0.0, 1.0, 5)

    with pytest.raises(ValueError, match="uniformly spaced"):
        solve_potential_dirichlet_2d(source, x, z)


def test_neumann_poisson_solver_recovers_cosine_mms_with_mean_gauge():
    x = np.linspace(0.0, 1.0, 49)
    z = np.linspace(0.0, 1.0, 51)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = np.cos(np.pi * x_grid) * np.cos(np.pi * z_grid)
    source = -2.0 * np.pi**2 * phi_exact

    phi = solve_potential_neumann_2d(source, x, z, mean_value=0.0)
    error = phi - phi_exact
    error -= np.mean(error)

    assert abs(np.mean(phi)) < 1.0e-12
    assert np.max(np.abs(error)) < 4.0e-4
    assert np.sqrt(np.mean(error**2)) < 2.0e-4


def test_neumann_poisson_solver_rejects_incompatible_source():
    x = np.linspace(0.0, 1.0, 5)
    z = np.linspace(0.0, 1.0, 5)
    source = np.ones((5, 5))

    with pytest.raises(ValueError, match="incompatible"):
        solve_potential_neumann_2d(source, x, z)


def test_neumann_poisson_solver_recovers_nonzero_gradient_polynomial_mms():
    x = np.linspace(0.0, 1.0, 33)
    z = np.linspace(0.0, 1.0, 35)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = (x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2
    phi_exact -= np.mean(phi_exact)
    source = np.full_like(phi_exact, 4.0)
    normal_gradient = {
        "x_min": np.ones_like(z),
        "x_max": np.ones_like(z),
        "z_min": np.ones_like(x),
        "z_max": np.ones_like(x),
    }

    phi = solve_potential_neumann_2d(source, x, z, mean_value=0.0, normal_gradient=normal_gradient)
    error = phi - phi_exact

    assert abs(np.mean(phi)) < 1.0e-12
    assert np.max(np.abs(error)) < 1.0e-10


def test_full_inductionless_reference_recovers_potential_and_current():
    x = np.linspace(0.0, 1.0, 49)
    z = np.linspace(0.0, 1.0, 51)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    sigma = 7.5
    b_y = 2.0
    magnetic_field = np.array([0.0, b_y, 0.0])

    phi_exact = (x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2
    phi_exact -= np.mean(phi_exact)
    streamfunction = 0.03 * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)

    emf = structured_gradient_2d(phi_exact, x, z) + target_current / sigma
    velocity = np.zeros_like(target_current)
    velocity[..., 0] = emf[..., 2] / b_y
    velocity[..., 2] = -emf[..., 0] / b_y
    np.testing.assert_allclose(motional_electric_field(velocity, magnetic_field), emf, rtol=1.0e-13, atol=1.0e-13)

    source = potential_source_from_motional_emf(emf, x, z)
    normal_gradient = insulating_wall_normal_gradient(emf, x, z)
    phi = solve_potential_neumann_2d(source, x, z, normal_gradient=normal_gradient, mean_value=0.0)
    current = reconstruct_inductionless_current_2d(phi, velocity, magnetic_field, sigma, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)

    np.testing.assert_allclose(phi, phi_exact, rtol=0.0, atol=1.0e-10)
    np.testing.assert_allclose(current, target_current, rtol=0.0, atol=1.0e-9)
    assert np.max(np.abs(residual[2:-2, 2:-2])) < 1.0e-9
    assert max(wall_current.values()) < 1.0e-9


def test_variable_conductivity_dirichlet_solver_recovers_smooth_mms():
    x = np.linspace(0.0, 1.0, 65)
    z = np.linspace(0.0, 1.0, 63)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    conductivity = 1.0 + 0.3 * x_grid + 0.2 * z_grid
    d_phi_dx = np.pi * np.cos(np.pi * x_grid) * np.sin(np.pi * z_grid)
    d_phi_dz = np.pi * np.sin(np.pi * x_grid) * np.cos(np.pi * z_grid)
    lap_phi = -2.0 * np.pi**2 * phi_exact
    source = conductivity * lap_phi + 0.3 * d_phi_dx + 0.2 * d_phi_dz

    phi = solve_variable_conductivity_dirichlet_2d(conductivity, source, x, z, boundary_values=0.0)
    error = phi - phi_exact

    assert np.max(np.abs(error)) < 3.0e-4
    assert np.sqrt(np.mean(error**2)) < 1.5e-4


def test_variable_conductivity_dirichlet_solver_rejects_nonpositive_sigma():
    x = np.linspace(0.0, 1.0, 5)
    z = np.linspace(0.0, 1.0, 5)
    source = np.zeros((5, 5))
    conductivity = np.ones((5, 5))
    conductivity[2, 2] = 0.0

    with pytest.raises(ValueError, match="conductivity"):
        solve_variable_conductivity_dirichlet_2d(conductivity, source, x, z)


def test_variable_conductivity_dirichlet_solver_handles_two_region_flux_continuity():
    x = np.linspace(0.0, 1.0, 64)
    z = np.linspace(0.0, 1.0, 17)
    x_grid, _ = np.meshgrid(x, z, indexing="xy")
    interface = 0.5
    sigma_left = 1.0
    sigma_right = 4.0
    flux = 1.0 / (interface / sigma_left + (1.0 - interface) / sigma_right)
    conductivity = np.where(x_grid < interface, sigma_left, sigma_right)
    phi_1d = np.where(
        x <= interface,
        flux * x / sigma_left,
        flux * interface / sigma_left + flux * (x - interface) / sigma_right,
    )
    phi_exact = np.broadcast_to(phi_1d, x_grid.shape)
    source = np.zeros_like(phi_exact)

    phi = solve_variable_conductivity_dirichlet_2d(conductivity, source, x, z, boundary_values=phi_exact)
    error = phi - phi_exact

    left_idx = x.size // 2 - 1
    right_idx = x.size // 2
    h = x[1] - x[0]
    sigma_face = 2.0 * sigma_left * sigma_right / (sigma_left + sigma_right)
    interface_flux = sigma_face * (phi[:, right_idx] - phi[:, left_idx]) / h

    assert np.max(np.abs(error)) < 1.0e-12
    np.testing.assert_allclose(interface_flux, flux, rtol=1.0e-12, atol=1.0e-12)


def test_variable_conductivity_neumann_solver_recovers_nonzero_flux_mms():
    x = np.linspace(0.0, 1.0, 49)
    z = np.linspace(0.0, 1.0, 51)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    phi_exact = (x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2
    phi_exact -= np.mean(phi_exact)
    conductivity = 1.0 + 0.3 * x_grid + 0.2 * z_grid
    source = 4.0 * conductivity + 0.6 * (x_grid - 0.5) + 0.4 * (z_grid - 0.5)
    normal_flux = {
        "x_min": conductivity[:, 0],
        "x_max": conductivity[:, -1],
        "z_min": conductivity[0, :],
        "z_max": conductivity[-1, :],
    }

    phi = solve_variable_conductivity_neumann_2d(
        conductivity,
        source,
        x,
        z,
        mean_value=0.0,
        normal_flux=normal_flux,
    )
    error = phi - phi_exact

    assert abs(np.mean(phi)) < 1.0e-12
    assert np.max(np.abs(error)) < 2.0e-4
    assert np.sqrt(np.mean(error**2)) < 1.0e-4


def test_variable_conductivity_neumann_solver_rejects_incompatible_flux():
    x = np.linspace(0.0, 1.0, 5)
    z = np.linspace(0.0, 1.0, 5)
    conductivity = np.ones((5, 5))
    source = np.ones((5, 5))

    with pytest.raises(ValueError, match="incompatible"):
        solve_variable_conductivity_neumann_2d(conductivity, source, x, z)


def test_full_variable_conductivity_inductionless_reference_is_charge_conserving():
    x = np.linspace(0.0, 1.0, 49)
    z = np.linspace(0.0, 1.0, 51)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    conductivity = 1.0 + 0.3 * x_grid + 0.2 * z_grid
    phi_exact = 0.02 * ((x_grid - 0.5) ** 2 + (z_grid - 0.5) ** 2)
    phi_exact -= np.mean(phi_exact)
    streamfunction = 0.03 * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)
    emf = structured_gradient_2d(phi_exact, x, z) + target_current / conductivity[..., None]
    magnetic_field = np.array([0.0, 2.0, 0.0])
    velocity = np.zeros_like(target_current)
    velocity[..., 0] = emf[..., 2] / magnetic_field[1]
    velocity[..., 2] = -emf[..., 0] / magnetic_field[1]

    motional = motional_electric_field(velocity, magnetic_field)
    source = weighted_potential_source_from_motional_emf(conductivity, motional, x, z)
    normal_flux = insulating_wall_normal_flux(conductivity, motional, x, z)
    phi = solve_variable_conductivity_neumann_2d(
        conductivity,
        source,
        x,
        z,
        mean_value=0.0,
        normal_flux=normal_flux,
    )
    current = reconstruct_inductionless_current_2d(phi, velocity, magnetic_field, conductivity, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)

    assert np.max(np.abs(phi - phi_exact)) < 5.0e-5
    assert np.max(np.abs(current - target_current)) < 5.0e-4
    assert np.max(np.abs(residual[2:-2, 2:-2])) < 2.0e-2
    assert max(wall_current.values()) < 5.0e-4


def test_two_region_full_inductionless_reference_checks_conservative_interface_current():
    nx = 64
    nz = 65
    x = np.linspace(0.0, 1.0, nx)
    z = np.linspace(0.0, 1.0, nz)
    x_grid, z_grid = np.meshgrid(x, z, indexing="xy")
    interface = 0.5
    sigma_left = 1.0
    sigma_right = 4.0
    flux = 0.7
    conductivity = np.where(x_grid < interface, sigma_left, sigma_right)
    phi_1d = np.where(
        x <= interface,
        flux * x / sigma_left,
        flux * interface / sigma_left + flux * (x - interface) / sigma_right,
    )
    phi_exact = np.broadcast_to(phi_1d, x_grid.shape).copy()
    phi_exact -= np.mean(phi_exact)
    streamfunction = 0.04 * np.sin(np.pi * x_grid) * np.sin(np.pi * z_grid)
    target_current = divergence_free_current_from_streamfunction(streamfunction, x, z)
    emf = np.zeros_like(target_current)
    emf[..., 0] = flux / conductivity + target_current[..., 0] / conductivity
    emf[..., 2] = target_current[..., 2] / conductivity
    magnetic_field = np.array([0.0, 3.0, 0.0])
    velocity = np.zeros_like(target_current)
    velocity[..., 0] = emf[..., 2] / magnetic_field[1]
    velocity[..., 2] = -emf[..., 0] / magnetic_field[1]
    motional = motional_electric_field(velocity, magnetic_field)

    source = weighted_potential_source_from_motional_emf(conductivity, motional, x, z)
    normal_flux = insulating_wall_normal_flux(conductivity, motional, x, z)
    phi = solve_variable_conductivity_neumann_2d(
        conductivity,
        source,
        x,
        z,
        mean_value=0.0,
        normal_flux=normal_flux,
    )
    current = reconstruct_inductionless_current_2d(phi, velocity, magnetic_field, conductivity, x, z)
    residual = charge_conservation_residual_2d(current, x, z)
    wall_current = wall_normal_current_extrema(current)

    left_idx = nx // 2 - 1
    right_idx = nx // 2
    bulk_mask = np.ones((nz, nx), dtype=bool)
    bulk_mask[:, left_idx - 1 : right_idx + 2] = False
    bulk_mask[:2, :] = False
    bulk_mask[-2:, :] = False
    bulk_mask[:, :2] = False
    bulk_mask[:, -2:] = False

    h = x[1] - x[0]
    sigma_face = 2.0 * sigma_left * sigma_right / (sigma_left + sigma_right)
    weighted_emf_x = conductivity * motional[..., 0]
    weighted_emf_face = 0.5 * (weighted_emf_x[:, left_idx] + weighted_emf_x[:, right_idx])
    interface_current = -sigma_face * (phi[:, right_idx] - phi[:, left_idx]) / h + weighted_emf_face
    target_interface_current = 0.5 * (target_current[:, left_idx, 0] + target_current[:, right_idx, 0])

    assert np.max(np.abs(phi - phi_exact)) < 1.0e-10
    assert np.max(np.abs((current - target_current)[bulk_mask])) < 1.0e-9
    assert np.max(np.abs(residual[bulk_mask])) < 1.0e-9
    assert max(wall_current.values()) < 1.0e-9
    np.testing.assert_allclose(interface_current, target_interface_current, rtol=0.0, atol=1.0e-9)
