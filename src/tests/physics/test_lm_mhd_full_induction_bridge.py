import importlib.util
from pathlib import Path

import numpy as np
import pytest


SCRIPT_PATH = (
    Path(__file__).resolve().parents[3] / "cases" / "lm_mhd_coupling" / "full_induction_low_rm_bridge.py"
)
SPEC = importlib.util.spec_from_file_location("full_induction_low_rm_bridge", SCRIPT_PATH)
bridge = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bridge)


def test_full_induction_bridge_recovers_inductionless_current():
    half_width = 0.05
    z = np.linspace(-half_width, half_width, 201)
    result = bridge.bridge_solution(
        z,
        velocity_scale=0.2,
        magnetic_field=5.0,
        conductivity=8.0e5,
        half_width=half_width,
    )

    assert result["rm"] == pytest.approx(bridge.MU0 * 8.0e5 * 0.2 * half_width)
    assert result["current_error"] < 3.0e-4
    assert abs(result["net_current"]) < 1.0e-9
    assert abs(result["wall_bx_mismatch"]) < 1.0e-14


def test_full_induction_bridge_induced_field_scales_with_rm():
    half_width = 0.05
    z = np.linspace(-half_width, half_width, 101)
    low = bridge.bridge_solution(z, 0.1, 5.0, 8.0e5, half_width)
    high = bridge.bridge_solution(z, 0.2, 5.0, 8.0e5, half_width)

    assert high["rm"] == pytest.approx(2.0 * low["rm"])
    assert high["induced_ratio"] == pytest.approx(2.0 * low["induced_ratio"], rel=1.0e-12)


def test_full_induction_bridge_transient_relaxes_to_steady_state():
    half_width = 0.05
    z = np.linspace(-half_width, half_width, 81)
    steady = bridge.bridge_solution(z, 0.2, 5.0, 8.0e5, half_width)

    _, error = bridge.transient_relaxation(
        z,
        steady["induced_bx"],
        velocity_scale=0.2,
        magnetic_field=5.0,
        conductivity=8.0e5,
        half_width=half_width,
        n_steps=80,
        t_end_factor=3.0,
    )

    assert error[-1] < 5.0e-4
    assert np.all(np.diff(error[1:]) <= 1.0e-12)
