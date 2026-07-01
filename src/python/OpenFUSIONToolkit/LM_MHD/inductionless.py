#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Inductionless LM-MHD reference utilities.

These helpers encode algebraic sign conventions and manufactured-solution
diagnostics for a future low-magnetic-Reynolds-number electric-potential mode.
They are not a solver and do not alter the Fortran MUG implementation.
'''
from __future__ import annotations

import numpy as np


def ohms_law_current_density(
    sigma: float | np.ndarray,
    velocity: np.ndarray,
    magnetic_field: np.ndarray,
    electric_potential_gradient: np.ndarray | None = None,
) -> np.ndarray:
    '''Return inductionless current density J = sigma(-grad(phi) + u x B).'''

    vel = _as_vector_array("velocity", velocity)
    bfield = _broadcast_vector("magnetic_field", magnetic_field, vel.shape)
    if electric_potential_gradient is None:
        grad_phi = np.zeros_like(vel)
    else:
        grad_phi = _broadcast_vector("electric_potential_gradient", electric_potential_gradient, vel.shape)
    sigma_arr = _positive_broadcast("sigma", sigma, vel.shape[:-1])
    return sigma_arr[..., None] * (-grad_phi + np.cross(vel, bfield))


def lorentz_force_density(current_density: np.ndarray, magnetic_field: np.ndarray) -> np.ndarray:
    '''Return Lorentz force density f = J x B.'''

    current = _as_vector_array("current_density", current_density)
    bfield = _broadcast_vector("magnetic_field", magnetic_field, current.shape)
    return np.cross(current, bfield)


def joule_heating_density(current_density: np.ndarray, sigma: float | np.ndarray) -> float | np.ndarray:
    '''Return q = |J|^2 / sigma, which is non-negative for positive sigma.'''

    current = _as_vector_array("current_density", current_density)
    sigma_arr = _positive_broadcast("sigma", sigma, current.shape[:-1])
    heating = np.sum(current * current, axis=-1) / sigma_arr
    return _maybe_scalar(heating)


def divergence_free_current_from_streamfunction(
    streamfunction: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> np.ndarray:
    '''Return a 2D divergence-free current from J_x = dpsi/dz, J_z = -dpsi/dx.

    `streamfunction` has shape `(nz, nx)` and `x`, `z` are one-dimensional
    coordinates. The returned array has shape `(nz, nx, 3)`, with the selected
    components populated and the remaining component set to zero.
    '''

    psi, xcoord, zcoord = _structured_scalar("streamfunction", streamfunction, x, z)
    comp_x, comp_z = _check_components(components)
    current = np.zeros(psi.shape + (3,), dtype=np.float64)
    current[..., comp_x] = np.gradient(psi, zcoord, axis=0, edge_order=2)
    current[..., comp_z] = -np.gradient(psi, xcoord, axis=1, edge_order=2)
    return current


def charge_conservation_residual_2d(
    current_density: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> np.ndarray:
    '''Return structured-grid residual div(J) = dJ_x/dx + dJ_z/dz.

    The current array must have shape `(nz, nx, 3)`. This diagnostic is useful
    for manufactured solutions and postprocessing; it is not a conservative
    finite-element divergence operator.
    '''

    current = _as_vector_array("current_density", current_density)
    if current.ndim != 3:
        raise ValueError("current_density must have shape (nz, nx, 3)")
    xcoord, zcoord = _structured_coordinates(x, z)
    if current.shape[:2] != (zcoord.size, xcoord.size):
        raise ValueError("current_density shape must match z and x coordinates")
    comp_x, comp_z = _check_components(components)
    d_jx_dx = np.gradient(current[..., comp_x], xcoord, axis=1, edge_order=2)
    d_jz_dz = np.gradient(current[..., comp_z], zcoord, axis=0, edge_order=2)
    return d_jx_dx + d_jz_dz


def wall_normal_current_extrema(
    current_density: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> dict[str, float]:
    '''Return max absolute normal current on rectangular grid boundaries.'''

    current = _as_vector_array("current_density", current_density)
    if current.ndim != 3:
        raise ValueError("current_density must have shape (nz, nx, 3)")
    comp_x, comp_z = _check_components(components)
    return {
        "x_min": float(np.max(np.abs(current[:, 0, comp_x]))),
        "x_max": float(np.max(np.abs(current[:, -1, comp_x]))),
        "z_min": float(np.max(np.abs(current[0, :, comp_z]))),
        "z_max": float(np.max(np.abs(current[-1, :, comp_z]))),
    }


def _as_vector_array(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape == () or arr.shape[-1] != 3:
        raise ValueError(f"{name} must have shape (..., 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _broadcast_vector(name: str, values: np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape == (3,):
        arr = np.broadcast_to(arr, target_shape)
    else:
        try:
            arr = np.broadcast_to(arr, target_shape)
        except ValueError as exc:
            raise ValueError(f"{name} must broadcast to {target_shape}") from exc
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _positive_broadcast(name: str, values: float | np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    try:
        arr = np.broadcast_to(arr, target_shape)
    except ValueError as exc:
        raise ValueError(f"{name} must broadcast to {target_shape}") from exc
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain positive finite values")
    return arr


def _structured_scalar(
    name: str,
    values: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"{name} must have shape (nz, nx)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    xcoord, zcoord = _structured_coordinates(x, z)
    if arr.shape != (zcoord.size, xcoord.size):
        raise ValueError(f"{name} shape must match z and x coordinates")
    return arr, xcoord, zcoord


def _structured_coordinates(x: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xcoord = np.asarray(x, dtype=np.float64)
    zcoord = np.asarray(z, dtype=np.float64)
    if xcoord.ndim != 1 or zcoord.ndim != 1:
        raise ValueError("x and z must be one-dimensional")
    if xcoord.size < 3 or zcoord.size < 3:
        raise ValueError("x and z must contain at least three points")
    if not np.all(np.isfinite(xcoord)) or not np.all(np.isfinite(zcoord)):
        raise ValueError("x and z must contain only finite values")
    if np.any(np.diff(xcoord) <= 0.0) or np.any(np.diff(zcoord) <= 0.0):
        raise ValueError("x and z must be strictly increasing")
    return xcoord, zcoord


def _check_components(components: tuple[int, int]) -> tuple[int, int]:
    if len(components) != 2:
        raise ValueError("components must contain two entries")
    comp_x, comp_z = components
    if comp_x == comp_z or comp_x not in (0, 1, 2) or comp_z not in (0, 1, 2):
        raise ValueError("components must be distinct vector component indices in 0..2")
    return comp_x, comp_z


def _maybe_scalar(values: np.ndarray) -> float | np.ndarray:
    if values.ndim == 0:
        return float(values)
    return values
