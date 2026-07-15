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


# --- Duct drag closures (LM-MHD gap-4 fix): validated vs resolved MUG runs -----
def _load_closures():
    import importlib.util
    from pathlib import Path
    p = Path(__file__).resolve().parents[2] / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "closures.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_closures_ducttest", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _u_nondim(alpha, fy=1e-2, nu=1.0, a=1.0):
    # reduced 0-D core: u = fy/alpha ; nondim u*nu/(fy a^2)
    return (fy / alpha) * nu / (fy * a**2)


def test_insulating_duct_core_scales_as_inverse_Ha():
    m = _load_closures()
    # insulating duct core ~ 1/Ha (Shercliff), NOT 1/Ha^2 (channel)
    for ha in (20.0, 100.0, 1000.0):
        assert _u_nondim(m.shercliff_duct_effective_drag(1.0, ha, 1.0)) == pytest.approx(1.0 / ha)
    # matches resolved MUG (c_wall=1e-4) to ~1%
    assert _u_nondim(m.shercliff_duct_effective_drag(1.0, 20.0, 1.0)) == pytest.approx(4.9831e-2, rel=0.01)
    assert _u_nondim(m.shercliff_duct_effective_drag(1.0, 100.0, 1.0)) == pytest.approx(9.9020e-3, rel=0.02)
    # the channel closure over-brakes an insulating duct by a factor of Ha
    ratio = m.hartmann_channel_effective_drag(1.0, 100.0, 1.0) / m.shercliff_duct_effective_drag(1.0, 100.0, 1.0)
    assert ratio == pytest.approx(100.0, rel=0.02)


def test_conducting_duct_core_matches_hunt_asymptote():
    m = _load_closures()
    # Hunt conducting-wall core -> (1 + 1/c)/Ha^2 at high Ha
    for c in (0.5, 1.0, 4.0):
        assert _u_nondim(m.hunt_duct_effective_drag(1.0, 1000.0, 1.0, c)) == pytest.approx((1.0 + 1.0 / c) / 1e6)
    # c=1 matches resolved MUG Hunt to <1% at Ha=100 (asymptote is high-Ha)
    assert _u_nondim(m.hunt_duct_effective_drag(1.0, 100.0, 1.0, 1.0)) == pytest.approx(1.995e-4, rel=0.01)
    # perfect-conductor limit c->inf recovers the channel Ha^2 drag
    perfect = m.hunt_duct_effective_drag(1.0, 100.0, 1.0, 1.0e8)
    channel = m.hartmann_channel_effective_drag(1.0, 100.0, 1.0)
    assert perfect == pytest.approx(channel, rel=0.02)


def test_duct_dispatcher_selects_wall_regime():
    m = _load_closures()
    ha = 100.0
    # conducting c=1 -> Hunt branch
    assert m.duct_effective_drag(1.0, ha, 1.0, 1.0) == pytest.approx(
        m.hunt_duct_effective_drag(1.0, ha, 1.0, 1.0))
    # insulating (c below the Ha^{-1/2} crossover) -> Shercliff branch
    assert m.duct_effective_drag(1.0, ha, 1.0, 1e-4) == pytest.approx(
        m.shercliff_duct_effective_drag(1.0, ha, 1.0))
    assert m.duct_effective_drag(1.0, ha, 1.0, 0.0) == pytest.approx(
        m.shercliff_duct_effective_drag(1.0, ha, 1.0))
