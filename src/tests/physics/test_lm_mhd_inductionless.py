import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        charge_conservation_residual_2d,
        divergence_free_current_from_streamfunction,
        insulating_wall_normal_gradient,
        joule_heating_density,
        lorentz_force_density,
        motional_electric_field,
        ohms_law_current_density,
        potential_source_from_motional_emf,
        reconstruct_inductionless_current_2d,
        solve_potential_dirichlet_2d,
        solve_potential_neumann_2d,
        structured_gradient_2d,
        wall_normal_current_extrema,
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
    charge_conservation_residual_2d = inductionless.charge_conservation_residual_2d
    divergence_free_current_from_streamfunction = inductionless.divergence_free_current_from_streamfunction
    insulating_wall_normal_gradient = inductionless.insulating_wall_normal_gradient
    joule_heating_density = inductionless.joule_heating_density
    lorentz_force_density = inductionless.lorentz_force_density
    motional_electric_field = inductionless.motional_electric_field
    ohms_law_current_density = inductionless.ohms_law_current_density
    potential_source_from_motional_emf = inductionless.potential_source_from_motional_emf
    reconstruct_inductionless_current_2d = inductionless.reconstruct_inductionless_current_2d
    solve_potential_dirichlet_2d = inductionless.solve_potential_dirichlet_2d
    solve_potential_neumann_2d = inductionless.solve_potential_neumann_2d
    structured_gradient_2d = inductionless.structured_gradient_2d
    wall_normal_current_extrema = inductionless.wall_normal_current_extrema


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
