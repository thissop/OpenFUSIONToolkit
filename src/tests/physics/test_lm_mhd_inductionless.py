import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        charge_conservation_residual_2d,
        divergence_free_current_from_streamfunction,
        joule_heating_density,
        lorentz_force_density,
        ohms_law_current_density,
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
    joule_heating_density = inductionless.joule_heating_density
    lorentz_force_density = inductionless.lorentz_force_density
    ohms_law_current_density = inductionless.ohms_law_current_density
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
