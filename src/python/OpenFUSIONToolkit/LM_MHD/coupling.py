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
from scipy import special

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


@dataclass(frozen=True)
class FeedbackMetrics:
    '''Reduced field-feedback metrics from axisymmetric toroidal currents.

    These metrics estimate the magnetic perturbation produced by currents in a
    `CouplingSnapshot` at user-chosen probe locations. They are triage numbers
    for deciding whether reverse coupling deserves a real TokaMaker solve, not
    replacements for that solve.
    '''

    total_toroidal_current: float
    absolute_toroidal_current: float
    current_centroid_rz: tuple[float, float]
    max_delta_b: float
    rms_delta_b: float
    max_feedback: float
    rms_feedback: float
    classification: str


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


def classify_feedback_ratio(
    ratio: float,
    feedback_threshold: float = 1.0e-3,
    strong_threshold: float = 1.0e-2,
) -> str:
    '''Classify reverse magnetic feedback from max |delta B|/|B_ref|.

    The thresholds are planning conventions. Crossing them is a prompt for a
    coupled sensitivity calculation, not a claim that the plasma equilibrium
    will move by the same fraction.
    '''

    _require_nonnegative("ratio", ratio)
    _require_positive("feedback_threshold", feedback_threshold)
    _require_positive("strong_threshold", strong_threshold)
    if strong_threshold <= feedback_threshold:
        raise ValueError("strong_threshold must exceed feedback_threshold")
    if ratio >= strong_threshold:
        return "strong-two-way-candidate"
    if ratio >= feedback_threshold:
        return "two-way-candidate"
    return "weak-feedback"


def axisymmetric_area_weights(snapshot: CouplingSnapshot) -> np.ndarray:
    '''Return poloidal-area weights from axisymmetric cell volumes.

    If `cell_volume` stores a toroidal volume for an axisymmetric cell, the
    corresponding poloidal cross-section area is dA = dV / (2 pi R). That dA is
    the weight needed to integrate toroidal current density, I = integral J_phi
    dA, for circular-loop field estimates.
    '''

    snap = snapshot.validated()
    if snap.cell_volume is None:
        raise ValueError("snapshot.cell_volume is required for axisymmetric area weights")
    radius = snap.points_rz[:, 0]
    if np.any(radius <= 0.0):
        raise ValueError("points_rz[:, 0] must be positive for axisymmetric area weights")
    return snap.cell_volume / (2.0 * np.pi * radius)


def integrated_toroidal_current(
    snapshot: CouplingSnapshot,
    area_weights: np.ndarray | None = None,
    absolute: bool = False,
) -> float:
    '''Return integral J_phi dA from a snapshot and poloidal-area weights.'''

    snap, weights = _snapshot_and_area_weights(snapshot, area_weights)
    current_density = snap.toroidal_current_density()
    if absolute:
        current_density = np.abs(current_density)
    return float(np.sum(current_density * weights))


def current_centroid_rz(
    snapshot: CouplingSnapshot,
    area_weights: np.ndarray | None = None,
    absolute: bool = True,
) -> tuple[float, float]:
    '''Return the J_phi-weighted current centroid in the poloidal plane.'''

    snap, weights = _snapshot_and_area_weights(snapshot, area_weights)
    current_weights = snap.toroidal_current_density() * weights
    if absolute:
        current_weights = np.abs(current_weights)
    denom = float(np.sum(current_weights))
    if abs(denom) == 0.0:
        raise ValueError("current centroid is undefined for zero weighted current")
    centroid = np.sum(snap.points_rz * current_weights[:, None], axis=0) / denom
    return float(centroid[0]), float(centroid[1])


def circular_loop_poloidal_field(
    loop_radius: float,
    loop_z: float,
    current: float,
    probe_rz: np.ndarray,
    mu: float = MU0,
    regularization: float = 0.0,
) -> np.ndarray:
    '''Return (B_R, B_Z) from an axisymmetric circular current loop.

    The expression is the standard Biot-Savart result for a circular loop,
    written with complete elliptic integrals. `regularization` optionally
    softens the loop/probe separation for diagnostic plots; it is not a
    substitute for a finite-element self-field treatment.
    '''

    _require_positive("loop_radius", loop_radius)
    _require_finite("loop_z", loop_z)
    _require_finite("current", current)
    _require_positive("mu", mu)
    _require_nonnegative("regularization", regularization)
    probes = _as_points_rz("probe_rz", probe_rz)
    radius = probes[:, 0]
    if np.any(radius < 0.0):
        raise ValueError("probe_rz[:, 0] must be non-negative")

    zeta = probes[:, 1] - loop_z
    zeta2 = zeta**2 + regularization**2
    separation2 = (loop_radius - radius) ** 2 + zeta2
    if regularization == 0.0 and np.any(separation2 == 0.0):
        raise ValueError("probe lies on the current loop; use regularization or move the probe")

    beta2 = (loop_radius + radius) ** 2 + zeta2
    beta = np.sqrt(beta2)
    m = np.zeros_like(radius)
    nonzero_beta = beta2 > 0.0
    m[nonzero_beta] = 4.0 * loop_radius * radius[nonzero_beta] / beta2[nonzero_beta]
    m = np.clip(m, 0.0, 1.0 - np.finfo(float).eps)

    ellip_k = special.ellipk(m)
    ellip_e = special.ellipe(m)
    common = mu * current / (2.0 * np.pi * beta)

    field_z = common * (
        ellip_k
        + (loop_radius**2 - radius**2 - zeta2) / separation2 * ellip_e
    )

    field_r = np.zeros_like(field_z)
    off_axis = radius > 0.0
    field_r[off_axis] = (
        common[off_axis]
        * zeta[off_axis]
        / radius[off_axis]
        * (
            -ellip_k[off_axis]
            + (loop_radius**2 + radius[off_axis] ** 2 + zeta2[off_axis])
            / separation2[off_axis]
            * ellip_e[off_axis]
        )
    )
    return np.column_stack((field_r, field_z))


def axisymmetric_toroidal_current_field(
    snapshot: CouplingSnapshot,
    probe_rz: np.ndarray,
    area_weights: np.ndarray | None = None,
    mu: float = MU0,
    regularization: float = 0.0,
) -> np.ndarray:
    '''Estimate (B_R, B_Z) at probes from sampled toroidal currents.

    Each sample contributes a circular loop with current dI = J_phi dA. This is
    a geometry-explicit feedback diagnostic for axisymmetric current snapshots.
    It is intentionally simple and should be replaced by the TokaMaker Green's
    function machinery for production coupled-equilibrium work.
    '''

    snap, weights = _snapshot_and_area_weights(snapshot, area_weights)
    probes = _as_points_rz("probe_rz", probe_rz)
    field = np.zeros((probes.shape[0], 2), dtype=np.float64)
    loop_currents = snap.toroidal_current_density() * weights
    for (radius, zval), loop_current in zip(snap.points_rz, loop_currents):
        if loop_current == 0.0:
            continue
        if radius <= 0.0:
            raise ValueError("current-loop radius must be positive")
        field += circular_loop_poloidal_field(
            radius,
            zval,
            float(loop_current),
            probes,
            mu=mu,
            regularization=regularization,
        )
    return field


def feedback_metrics_from_snapshot(
    snapshot: CouplingSnapshot,
    probe_rz: np.ndarray,
    reference_field: float,
    area_weights: np.ndarray | None = None,
    regularization: float = 0.0,
) -> FeedbackMetrics:
    '''Return reduced magnetic-feedback metrics for a current snapshot.'''

    _require_positive("reference_field", abs(reference_field))
    field = axisymmetric_toroidal_current_field(
        snapshot,
        probe_rz,
        area_weights=area_weights,
        regularization=regularization,
    )
    field_mag = np.linalg.norm(field, axis=1)
    max_delta_b = float(np.max(field_mag))
    rms_delta_b = float(np.sqrt(np.mean(field_mag**2)))
    max_ratio = float(max_delta_b / abs(reference_field))
    rms_ratio = float(rms_delta_b / abs(reference_field))
    return FeedbackMetrics(
        total_toroidal_current=integrated_toroidal_current(snapshot, area_weights=area_weights),
        absolute_toroidal_current=integrated_toroidal_current(snapshot, area_weights=area_weights, absolute=True),
        current_centroid_rz=current_centroid_rz(snapshot, area_weights=area_weights, absolute=True),
        max_delta_b=max_delta_b,
        rms_delta_b=rms_delta_b,
        max_feedback=max_ratio,
        rms_feedback=rms_ratio,
        classification=classify_feedback_ratio(max_ratio),
    )


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


def _as_points_rz(name: str, values: Any) -> np.ndarray:
    arr = _as_float_array(name, values, 2)
    if arr.shape[1] != 2:
        raise ValueError(f"{name} must have shape (n, 2)")
    return arr


def _area_weights(values: Any, size: int) -> np.ndarray:
    weights = _as_float_array("area_weights", values, 1)
    if weights.shape != (size,):
        raise ValueError("area_weights must have shape (n,)")
    if np.any(weights < 0.0):
        raise ValueError("area_weights entries must be non-negative")
    return weights


def _snapshot_and_area_weights(
    snapshot: CouplingSnapshot,
    area_weights: np.ndarray | None,
) -> tuple[CouplingSnapshot, np.ndarray]:
    snap = snapshot.validated()
    if area_weights is None:
        weights = axisymmetric_area_weights(snap)
    else:
        weights = _area_weights(area_weights, snap.points_rz.shape[0])
    return snap, weights


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
