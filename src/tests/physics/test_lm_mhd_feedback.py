import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        MU0,
        CouplingSnapshot,
        axisymmetric_area_weights,
        axisymmetric_toroidal_current_field,
        circular_loop_poloidal_field,
        classify_feedback_ratio,
        current_centroid_rz,
        feedback_metrics_from_snapshot,
        integrated_toroidal_current,
    )
except FileNotFoundError:
    # Source-tree fallback for Mac-side development shells where the compiled
    # OFT shared library is not present next to the Python package.
    import importlib.util
    import sys
    from pathlib import Path

    coupling_path = Path(__file__).resolve().parents[2] / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = coupling
    spec.loader.exec_module(coupling)
    MU0 = coupling.MU0
    CouplingSnapshot = coupling.CouplingSnapshot
    axisymmetric_area_weights = coupling.axisymmetric_area_weights
    axisymmetric_toroidal_current_field = coupling.axisymmetric_toroidal_current_field
    circular_loop_poloidal_field = coupling.circular_loop_poloidal_field
    classify_feedback_ratio = coupling.classify_feedback_ratio
    current_centroid_rz = coupling.current_centroid_rz
    feedback_metrics_from_snapshot = coupling.feedback_metrics_from_snapshot
    integrated_toroidal_current = coupling.integrated_toroidal_current


def _snapshot(j_phi, points_rz=None, cell_volume=None):
    if points_rz is None:
        points_rz = np.array([[1.0, -0.1], [1.5, 0.2]])
    current_density = np.zeros((len(points_rz), 3))
    current_density[:, 1] = np.asarray(j_phi, dtype=float)
    return CouplingSnapshot(
        points_rz=np.asarray(points_rz, dtype=float),
        current_density=current_density,
        time=0.0,
        cell_volume=cell_volume,
    )


def test_axisymmetric_area_weights_convert_volume_to_poloidal_area():
    points = np.array([[1.0, 0.0], [1.4, 0.2], [1.8, -0.2]])
    area = np.array([1.0e-4, 2.0e-4, 3.0e-4])
    volume = 2.0 * np.pi * points[:, 0] * area
    snap = _snapshot([1.0, 2.0, 3.0], points_rz=points, cell_volume=volume)

    np.testing.assert_allclose(axisymmetric_area_weights(snap), area)


def test_integrated_current_and_centroid_use_explicit_area_weights():
    points = np.array([[1.0, -0.2], [2.0, 0.3]])
    area = np.array([1.0e-4, 3.0e-4])
    snap = _snapshot([2.0e5, -1.0e5], points_rz=points)

    assert integrated_toroidal_current(snap, area_weights=area) == pytest.approx(-10.0)
    assert integrated_toroidal_current(snap, area_weights=area, absolute=True) == pytest.approx(50.0)
    centroid = current_centroid_rz(snap, area_weights=area, absolute=True)
    expected = (20.0 * points[0] + 30.0 * points[1]) / 50.0
    np.testing.assert_allclose(centroid, expected)


def test_circular_loop_field_matches_axis_closed_form():
    loop_radius = 1.25
    loop_z = -0.15
    current = 17.0
    probes = np.array([[0.0, -0.15], [0.0, 0.2], [0.0, -0.7]])

    field = circular_loop_poloidal_field(loop_radius, loop_z, current, probes)
    zeta = probes[:, 1] - loop_z
    expected_bz = MU0 * current * loop_radius**2 / (2.0 * (loop_radius**2 + zeta**2) ** 1.5)

    np.testing.assert_allclose(field[:, 0], 0.0, atol=1.0e-18)
    np.testing.assert_allclose(field[:, 1], expected_bz, rtol=1.0e-12, atol=0.0)


def test_axisymmetric_toroidal_current_field_scales_linearly():
    points = np.array([[1.0, 0.0]])
    area = np.array([2.0e-4])
    probe = np.array([[0.7, 0.3], [1.6, -0.25]])
    snap = _snapshot([1.5e5], points_rz=points)
    snap_double = _snapshot([3.0e5], points_rz=points)

    field = axisymmetric_toroidal_current_field(snap, probe, area_weights=area)
    field_double = axisymmetric_toroidal_current_field(snap_double, probe, area_weights=area)

    np.testing.assert_allclose(field_double, 2.0 * field, rtol=1.0e-12, atol=0.0)


def test_feedback_metrics_are_finite_and_classified():
    points = np.array([[1.0, -0.1], [1.2, 0.1]])
    area = np.array([1.5e-4, 1.5e-4])
    snap = _snapshot([2.0e5, 1.0e5], points_rz=points)
    probes = np.array([[0.8, 0.0], [1.6, 0.0], [2.0, 0.3]])

    metrics = feedback_metrics_from_snapshot(snap, probes, reference_field=5.0, area_weights=area)

    assert metrics.total_toroidal_current == pytest.approx(45.0)
    assert metrics.absolute_toroidal_current == pytest.approx(45.0)
    assert np.isfinite(metrics.max_delta_b)
    assert np.isfinite(metrics.rms_delta_b)
    assert metrics.max_delta_b >= metrics.rms_delta_b > 0.0
    assert metrics.max_feedback == pytest.approx(metrics.max_delta_b / 5.0)
    assert metrics.classification == classify_feedback_ratio(metrics.max_feedback)


def test_feedback_classifier_thresholds_are_ordered():
    assert classify_feedback_ratio(5.0e-4) == "weak-feedback"
    assert classify_feedback_ratio(2.0e-3) == "two-way-candidate"
    assert classify_feedback_ratio(2.0e-2) == "strong-two-way-candidate"


def test_loop_singularity_requires_regularization():
    probe_on_loop = np.array([[1.0, 0.0]])
    with pytest.raises(ValueError, match="probe lies on the current loop"):
        circular_loop_poloidal_field(1.0, 0.0, 1.0, probe_on_loop)
    field = circular_loop_poloidal_field(1.0, 0.0, 1.0, probe_on_loop, regularization=1.0e-3)
    assert np.all(np.isfinite(field))
