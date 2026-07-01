#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Liquid-metal MHD regime metrics and a minimal coupling snapshot schema.

The functions in this module are deliberately small and formula-level. They are
intended for pre-solve triage, order-of-magnitude red-team checks, and I/O
contract tests, not as a substitute for a resolved or inductionless MHD solve.
'''
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import h5py
import numpy as np

MU0 = 4.0e-7 * np.pi
SCHEMA_VERSION = "lm_mhd_coupling_snapshot_v0"


@dataclass(frozen=True)
class LiquidMetalProperties:
    '''Uniform liquid-metal material properties used in MHD scale estimates.

    Parameters are SI: density in kg/m^3, kinematic viscosity in m^2/s, and
    electrical conductivity in S/m. The bundled representative values are only
    starting points for regime scans; design calculations should pass a
    temperature-specific property set.
    '''

    name: str
    density: float
    kinematic_viscosity: float
    electrical_conductivity: float
    temperature_c: float | None = None
    note: str = ""

    def __post_init__(self) -> None:
        _require_positive("density", self.density)
        _require_positive("kinematic_viscosity", self.kinematic_viscosity)
        _require_positive("electrical_conductivity", self.electrical_conductivity)


@dataclass(frozen=True)
class CouplingSnapshot:
    '''Mesh-sampled current data intended for MUG-to-TokaMaker feedback work.

    `points_rz` are sample locations in the poloidal plane, with shape `(n, 2)`.
    `current_density` has shape `(n, 3)` and components `(R, phi, Z)` in A/m^2.
    `time` is seconds. Optional `region_id` and `cell_volume` arrays are
    sample-wise metadata for projection/integration in later coupling work.
    '''

    points_rz: np.ndarray
    current_density: np.ndarray
    time: float
    region_id: np.ndarray | None = None
    cell_volume: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def validated(self) -> "CouplingSnapshot":
        points = _as_float_array("points_rz", self.points_rz, 2)
        currents = _as_float_array("current_density", self.current_density, 2)
        if points.shape[1] != 2:
            raise ValueError("points_rz must have shape (n, 2)")
        if currents.shape != (points.shape[0], 3):
            raise ValueError("current_density must have shape (n, 3)")
        _require_finite("time", self.time)
        if self.region_id is None:
            region_id = None
        else:
            region_id = np.asarray(self.region_id, dtype=np.int32)
            if region_id.shape != (points.shape[0],):
                raise ValueError("region_id must have shape (n,)")
        if self.cell_volume is None:
            cell_volume = None
        else:
            cell_volume = _as_float_array("cell_volume", self.cell_volume, 1)
            if cell_volume.shape != (points.shape[0],):
                raise ValueError("cell_volume must have shape (n,)")
            if np.any(cell_volume <= 0.0):
                raise ValueError("cell_volume entries must be positive")
        return CouplingSnapshot(
            points_rz=points,
            current_density=currents,
            time=float(self.time),
            region_id=region_id,
            cell_volume=cell_volume,
            metadata=dict(self.metadata),
        )

    def toroidal_current_density(self) -> np.ndarray:
        '''Return the toroidal (`phi`) current component in A/m^2.'''

        snap = self.validated()
        return snap.current_density[:, 1].copy()


def representative_liquid_metals() -> dict[str, LiquidMetalProperties]:
    '''Return representative PbLi and Li properties for exploratory scans.

    The values are intentionally rounded because published correlations vary
    with alloy composition and temperature. They are appropriate for sensitivity
    studies and dimensional checks only.
    '''

    return {
        "PbLi": LiquidMetalProperties(
            name="PbLi, representative near 500 C",
            density=9.5e3,
            kinematic_viscosity=1.8e-7,
            electrical_conductivity=0.8e6,
            temperature_c=500.0,
            note="Rounded fusion-blanket PbLi scale; replace with a chosen property correlation before design use.",
        ),
        "Li": LiquidMetalProperties(
            name="Li, representative near 500 C",
            density=5.0e2,
            kinematic_viscosity=9.0e-7,
            electrical_conductivity=3.2e6,
            temperature_c=500.0,
            note="Rounded liquid-lithium scale; replace with a chosen property correlation before design use.",
        ),
    }


def reynolds_number(velocity: float | np.ndarray, length: float, nu: float) -> float | np.ndarray:
    '''Return Re = U L / nu.'''

    _require_positive("length", length)
    _require_positive("nu", nu)
    return np.asarray(velocity) * length / nu


def magnetic_reynolds_number(
    velocity: float | np.ndarray,
    length: float,
    sigma: float,
    mu: float = MU0,
) -> float | np.ndarray:
    '''Return Rm = mu sigma U L.'''

    _require_positive("length", length)
    _require_positive("sigma", sigma)
    _require_positive("mu", mu)
    return mu * sigma * np.asarray(velocity) * length


def hartmann_number(
    magnetic_field: float | np.ndarray,
    length: float,
    sigma: float,
    density: float,
    nu: float,
) -> float | np.ndarray:
    '''Return Ha = B L sqrt(sigma / (rho nu)).'''

    _require_positive("length", length)
    _require_positive("sigma", sigma)
    _require_positive("density", density)
    _require_positive("nu", nu)
    return np.asarray(magnetic_field) * length * np.sqrt(sigma / (density * nu))


def interaction_parameter(
    magnetic_field: float | np.ndarray,
    length: float,
    velocity: float | np.ndarray,
    sigma: float,
    density: float,
) -> float | np.ndarray:
    '''Return N = sigma B^2 L / (rho U), the Lorentz-to-inertia scale.'''

    _require_positive("length", length)
    _require_positive("sigma", sigma)
    _require_positive("density", density)
    velocity_arr = np.asarray(velocity)
    if np.any(velocity_arr <= 0.0):
        raise ValueError("velocity must be positive for interaction_parameter")
    return sigma * np.asarray(magnetic_field) ** 2 * length / (density * velocity_arr)


def magnetic_diffusion_time(length: float, sigma: float, mu: float = MU0) -> float:
    '''Return tau_m = mu sigma L^2.'''

    _require_positive("length", length)
    _require_positive("sigma", sigma)
    _require_positive("mu", mu)
    return float(mu * sigma * length**2)


def flow_time(velocity: float, length: float) -> float:
    '''Return tau_u = L / U.'''

    _require_positive("velocity", velocity)
    _require_positive("length", length)
    return float(length / velocity)


def alfven_speed(magnetic_field: float, density: float, mu: float = MU0) -> float:
    '''Return v_A = B / sqrt(mu rho).'''

    _require_positive("density", density)
    _require_positive("mu", mu)
    return float(abs(magnetic_field) / np.sqrt(mu * density))


def alfven_mach_number(velocity: float, magnetic_field: float, density: float, mu: float = MU0) -> float:
    '''Return M_A = U / v_A.'''

    _require_positive("velocity", velocity)
    return float(velocity / alfven_speed(magnetic_field, density, mu))


def motional_current_density(
    sigma: float,
    velocity: float | np.ndarray,
    magnetic_field: float | np.ndarray,
    closure_factor: float = 1.0,
) -> float | np.ndarray:
    '''Return the scale j ~ closure_factor * sigma U B.

    `closure_factor=1` is a conservative upper-bound scale. In real ducts the
    electric potential and wall/current closure can strongly reduce or redirect
    this current, so this is a triage number rather than a prediction.
    '''

    _require_positive("sigma", sigma)
    if closure_factor < 0.0:
        raise ValueError("closure_factor must be non-negative")
    return closure_factor * sigma * np.asarray(velocity) * np.asarray(magnetic_field)


def sheet_delta_b(
    current_density: float | np.ndarray,
    conducting_thickness: float,
    mu: float = MU0,
) -> float | np.ndarray:
    '''Return the sheet-current field scale delta B ~ mu K / 2, K = j t.'''

    _require_positive("conducting_thickness", conducting_thickness)
    _require_positive("mu", mu)
    return 0.5 * mu * np.asarray(current_density) * conducting_thickness


def feedback_ratio(delta_b: float | np.ndarray, reference_field: float) -> float | np.ndarray:
    '''Return |delta B| / |B_ref|.'''

    _require_positive("reference_field", abs(reference_field))
    return np.abs(np.asarray(delta_b)) / abs(reference_field)


def blanket_feedback_ratio_sheet(
    props: LiquidMetalProperties,
    velocity: float | np.ndarray,
    magnetic_field: float | np.ndarray,
    conducting_thickness: float,
    reference_field: float,
    closure_factor: float = 1.0,
) -> float | np.ndarray:
    '''Upper-bound two-way feedback scale from a motional sheet current.'''

    j = motional_current_density(
        props.electrical_conductivity,
        velocity,
        magnetic_field,
        closure_factor=closure_factor,
    )
    return feedback_ratio(sheet_delta_b(j, conducting_thickness), reference_field)


def classify_coupling(
    *,
    rm: float,
    feedback: float,
    interaction: float,
    rm_threshold: float = 0.1,
    feedback_threshold: float = 1.0e-3,
    interaction_threshold: float = 1.0,
) -> str:
    '''Classify whether one-way, two-way, and drag physics deserve attention.

    Returns one of:
    - "low": low-Rm one-way field is likely adequate for first-pass scans.
    - "drag-important": Lorentz forces are strong but magnetic feedback is small.
    - "two-way-candidate": blanket currents may perturb the applied field enough
      to warrant a coupled calculation.
    - "full-induction-candidate": Rm is high enough that induced-field evolution
      should be questioned rather than assumed inductionless.
    '''

    _require_nonnegative("rm", rm)
    _require_nonnegative("feedback", feedback)
    _require_nonnegative("interaction", interaction)
    if rm >= rm_threshold:
        return "full-induction-candidate"
    if feedback >= feedback_threshold:
        return "two-way-candidate"
    if interaction >= interaction_threshold:
        return "drag-important"
    return "low"


def write_coupling_snapshot(filename: str, snapshot: CouplingSnapshot) -> None:
    '''Write a validated current-feedback snapshot to HDF5.'''

    snap = snapshot.validated()
    with h5py.File(filename, "w") as h5:
        root = h5.create_group("lm_mhd")
        root.attrs["schema_version"] = SCHEMA_VERSION
        root.attrs["time_units"] = "s"
        root.attrs["coordinate_system"] = "axisymmetric_R_phi_Z"
        root.attrs["current_density_units"] = "A m^-2"
        for key, value in snap.metadata.items():
            if _hdf5_attr_supported(value):
                root.attrs[key] = value
        root.create_dataset("time", data=snap.time)
        root.create_dataset("points_RZ", data=snap.points_rz)
        root.create_dataset("current_density_R_phi_Z", data=snap.current_density)
        if snap.region_id is not None:
            root.create_dataset("region_id", data=snap.region_id)
        if snap.cell_volume is not None:
            dset = root.create_dataset("cell_volume", data=snap.cell_volume)
            dset.attrs["units"] = "m^3"


def read_coupling_snapshot(filename: str) -> CouplingSnapshot:
    '''Read and validate an HDF5 current-feedback snapshot.'''

    with h5py.File(filename, "r") as h5:
        if "lm_mhd" not in h5:
            raise ValueError("missing /lm_mhd group")
        root = h5["lm_mhd"]
        schema = root.attrs.get("schema_version", "")
        if _decode_attr(schema) != SCHEMA_VERSION:
            raise ValueError(f"unsupported schema_version {schema!r}")
        required = ("time", "points_RZ", "current_density_R_phi_Z")
        for name in required:
            if name not in root:
                raise ValueError(f"missing /lm_mhd/{name}")
        metadata = {key: _decode_attr(value) for key, value in root.attrs.items()}
        snap = CouplingSnapshot(
            points_rz=np.asarray(root["points_RZ"]),
            current_density=np.asarray(root["current_density_R_phi_Z"]),
            time=float(np.asarray(root["time"])),
            region_id=np.asarray(root["region_id"]) if "region_id" in root else None,
            cell_volume=np.asarray(root["cell_volume"]) if "cell_volume" in root else None,
            metadata=metadata,
        )
    return snap.validated()


def _as_float_array(name: str, values: Any, ndim: int) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _require_positive(name: str, value: float) -> None:
    _require_finite(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")


def _require_nonnegative(name: str, value: float) -> None:
    _require_finite(name, value)
    if value < 0.0:
        raise ValueError(f"{name} must be non-negative")


def _require_finite(name: str, value: float) -> None:
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")


def _hdf5_attr_supported(value: Any) -> bool:
    return isinstance(value, (str, bytes, int, float, np.integer, np.floating, np.bool_))


def _decode_attr(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value
