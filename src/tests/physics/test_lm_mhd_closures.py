import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        cells_per_hartmann_layer,
        classify_hartmann_resolution,
        drag_power_density,
        hartmann_channel_effective_drag,
        hartmann_layer_thickness,
        hartmann_mean_over_centerline,
        implicit_drag_update,
        parallel_velocity,
        perpendicular_velocity,
    )
except FileNotFoundError:
    import importlib.util
    import sys
    from pathlib import Path

    closures_path = Path(__file__).resolve().parents[2] / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "closures.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_closures", closures_path)
    closures = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = closures
    spec.loader.exec_module(closures)
    cells_per_hartmann_layer = closures.cells_per_hartmann_layer
    classify_hartmann_resolution = closures.classify_hartmann_resolution
    drag_power_density = closures.drag_power_density
    hartmann_channel_effective_drag = closures.hartmann_channel_effective_drag
    hartmann_layer_thickness = closures.hartmann_layer_thickness
    hartmann_mean_over_centerline = closures.hartmann_mean_over_centerline
    implicit_drag_update = closures.implicit_drag_update
    parallel_velocity = closures.parallel_velocity
    perpendicular_velocity = closures.perpendicular_velocity


def test_hartmann_layer_thickness_and_cells_per_layer():
    assert hartmann_layer_thickness(0.2, 20.0) == pytest.approx(0.01)
    assert cells_per_hartmann_layer(0.0025, 0.2, 20.0) == pytest.approx(4.0)


def test_resolution_classifier_thresholds():
    a = 1.0
    ha = 10.0
    assert classify_hartmann_resolution(a / 80.0, a, ha) == "resolved"
    assert classify_hartmann_resolution(a / 50.0, a, ha) == "marginal"
    assert classify_hartmann_resolution(a / 20.0, a, ha) == "underresolved"


def test_effective_drag_limits_and_positive_values():
    a = 0.5
    nu = 2.0e-7
    assert hartmann_channel_effective_drag(a, 0.0, nu) == pytest.approx(3.0 * nu / a**2)
    assert hartmann_channel_effective_drag(a, 1.0e-7, nu) == pytest.approx(3.0 * nu / a**2, rel=1.0e-12)

    ha = np.array([1.0, 5.0, 50.0])
    alpha = hartmann_channel_effective_drag(a, ha, nu)
    assert np.all(np.isfinite(alpha))
    assert np.all(alpha > 0.0)
    assert alpha[-1] == pytest.approx(nu * ha[-1] ** 2 / a**2, rel=0.03)
    assert hartmann_mean_over_centerline(0.0) == pytest.approx(2.0 / 3.0)


def test_drag_update_preserves_parallel_velocity_and_damps_perpendicular():
    velocities = np.array([[1.0, 2.0, 3.0], [-2.0, 0.5, 1.5]])
    bhat = np.array([0.0, 2.0, 0.0])
    updated = implicit_drag_update(velocities, alpha=4.0, dt=0.25, bhat=bhat)
    np.testing.assert_allclose(parallel_velocity(updated, bhat), parallel_velocity(velocities, bhat))
    np.testing.assert_allclose(
        perpendicular_velocity(updated, bhat),
        perpendicular_velocity(velocities, bhat) / 2.0,
    )


def test_drag_power_is_non_positive_and_parallel_flow_has_zero_power():
    velocities = np.array([[1.0, 2.0, 3.0], [0.0, -2.0, 0.0], [4.0, 0.0, 0.0]])
    bhat = np.array([0.0, 1.0, 0.0])
    power = drag_power_density(velocities, alpha=3.0, bhat=bhat, density=7.0)
    assert np.all(power <= 0.0)
    assert power[1] == pytest.approx(0.0)
    assert power[0] < 0.0
    assert power[2] < 0.0


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError, match="hartmann"):
        hartmann_layer_thickness(1.0, 0.0)
    with pytest.raises(ValueError, match="bhat"):
        implicit_drag_update(np.ones(3), alpha=1.0, dt=0.1, bhat=np.zeros(3))
