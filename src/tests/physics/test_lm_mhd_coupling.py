import numpy as np
import pytest

try:
    from OpenFUSIONToolkit.LM_MHD import (
        CouplingSnapshot,
        MU0,
        blanket_feedback_ratio_sheet,
        classify_coupling,
        flow_time,
        hartmann_number,
        interaction_parameter,
        magnetic_diffusion_time,
        magnetic_reynolds_number,
        read_coupling_snapshot,
        representative_liquid_metals,
        reynolds_number,
        sheet_delta_b,
        structured_snapshot_current_grid,
        write_coupling_snapshot,
    )
except FileNotFoundError:
    # Source-tree fallback for Mac-side development shells where the compiled
    # OFT shared library is not present next to the Python package.
    import importlib.util
    from pathlib import Path

    coupling_path = Path(__file__).resolve().parents[2] / "python" / "OpenFUSIONToolkit" / "LM_MHD" / "coupling.py"
    spec = importlib.util.spec_from_file_location("lm_mhd_coupling", coupling_path)
    coupling = importlib.util.module_from_spec(spec)
    import sys

    sys.modules[spec.name] = coupling
    spec.loader.exec_module(coupling)
    CouplingSnapshot = coupling.CouplingSnapshot
    MU0 = coupling.MU0
    blanket_feedback_ratio_sheet = coupling.blanket_feedback_ratio_sheet
    classify_coupling = coupling.classify_coupling
    flow_time = coupling.flow_time
    hartmann_number = coupling.hartmann_number
    interaction_parameter = coupling.interaction_parameter
    magnetic_diffusion_time = coupling.magnetic_diffusion_time
    magnetic_reynolds_number = coupling.magnetic_reynolds_number
    read_coupling_snapshot = coupling.read_coupling_snapshot
    representative_liquid_metals = coupling.representative_liquid_metals
    reynolds_number = coupling.reynolds_number
    sheet_delta_b = coupling.sheet_delta_b
    structured_snapshot_current_grid = coupling.structured_snapshot_current_grid
    write_coupling_snapshot = coupling.write_coupling_snapshot


def test_dimensionless_identities():
    props = representative_liquid_metals()["PbLi"]
    B = 5.0
    U = 0.2
    L = 0.05

    ha = hartmann_number(B, L, props.electrical_conductivity, props.density, props.kinematic_viscosity)
    re = reynolds_number(U, L, props.kinematic_viscosity)
    n_interaction = interaction_parameter(B, L, U, props.electrical_conductivity, props.density)
    rm = magnetic_reynolds_number(U, L, props.electrical_conductivity)

    assert n_interaction == pytest.approx(ha**2 / re)
    assert rm == pytest.approx(magnetic_diffusion_time(L, props.electrical_conductivity) / flow_time(U, L))


def test_sheet_feedback_scale_is_motional_upper_bound():
    props = representative_liquid_metals()["Li"]
    ratio = blanket_feedback_ratio_sheet(
        props,
        velocity=0.1,
        magnetic_field=5.0,
        conducting_thickness=0.02,
        reference_field=5.0,
        closure_factor=0.25,
    )
    expected_j = 0.25 * props.electrical_conductivity * 0.1 * 5.0
    expected = abs(0.5 * MU0 * expected_j * 0.02) / 5.0
    assert ratio == pytest.approx(expected)


def test_classification_thresholds_are_ordered():
    assert classify_coupling(rm=0.01, feedback=1e-5, interaction=0.1) == "low"
    assert classify_coupling(rm=0.01, feedback=1e-5, interaction=10.0) == "drag-important"
    assert classify_coupling(rm=0.01, feedback=2e-3, interaction=0.1) == "two-way-candidate"
    assert classify_coupling(rm=0.2, feedback=1e-5, interaction=0.1) == "full-induction-candidate"


def test_hdf5_snapshot_roundtrip(tmp_path):
    snap = CouplingSnapshot(
        points_rz=np.array([[1.0, -0.2], [1.1, 0.0], [1.2, 0.2]]),
        current_density=np.array([[0.0, 1.0e5, 0.0], [1.0, 2.0e5, 3.0], [0.0, -1.0e5, 0.0]]),
        time=1.25e-3,
        region_id=np.array([3, 3, 4]),
        cell_volume=np.array([1.0e-4, 1.2e-4, 1.0e-4]),
        metadata={"source": "unit-test"},
    )
    path = tmp_path / "snapshot.h5"
    write_coupling_snapshot(str(path), snap)
    out = read_coupling_snapshot(str(path))

    np.testing.assert_allclose(out.points_rz, snap.points_rz)
    np.testing.assert_allclose(out.current_density, snap.current_density)
    np.testing.assert_allclose(out.toroidal_current_density(), snap.current_density[:, 1])
    assert out.time == pytest.approx(snap.time)
    assert out.metadata["schema_version"] == "lm_mhd_coupling_snapshot_v0"


def test_structured_snapshot_current_grid_sorts_shuffled_points():
    r = np.array([1.0, 1.5, 2.0])
    z = np.array([-0.25, 0.25])
    r_grid, z_grid = np.meshgrid(r, z, indexing="xy")
    current = np.zeros(r_grid.shape + (3,))
    current[..., 0] = r_grid
    current[..., 1] = 2.0 * z_grid
    current[..., 2] = r_grid + z_grid
    order = np.array([4, 0, 5, 2, 1, 3])
    snap = CouplingSnapshot(
        points_rz=np.column_stack((r_grid.ravel(), z_grid.ravel()))[order],
        current_density=current.reshape(-1, 3)[order],
        time=0.0,
    )

    out_r, out_z, out_current = structured_snapshot_current_grid(snap)

    np.testing.assert_allclose(out_r, r)
    np.testing.assert_allclose(out_z, z)
    np.testing.assert_allclose(out_current, current)


def test_structured_snapshot_current_grid_rejects_incomplete_grid():
    snap = CouplingSnapshot(
        points_rz=np.array([[1.0, 0.0], [1.5, 0.0], [1.0, 0.5]]),
        current_density=np.zeros((3, 3)),
        time=0.0,
    )
    with pytest.raises(ValueError, match="structured"):
        structured_snapshot_current_grid(snap)


def test_snapshot_rejects_inconsistent_shapes():
    snap = CouplingSnapshot(
        points_rz=np.zeros((2, 2)),
        current_density=np.zeros((3, 3)),
        time=0.0,
    )
    with pytest.raises(ValueError, match="current_density"):
        snap.validated()


def test_sheet_delta_b_rejects_nonphysical_thickness():
    with pytest.raises(ValueError, match="conducting_thickness"):
        sheet_delta_b(1.0e5, 0.0)
