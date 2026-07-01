#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''Reduced LM-MHD closure helpers for unresolved Hartmann-layer diagnostics.

These functions are analytic or algebraic checks for planning/calibration. They
do not replace a resolved MHD solve and should not be used to claim pointwise
Hartmann-layer agreement.
'''
from __future__ import annotations

import numpy as np


def hartmann_layer_thickness(half_width: float, hartmann: float | np.ndarray) -> float | np.ndarray:
    '''Return the Hartmann layer scale delta_H = a / Ha.'''

    _require_positive("half_width", half_width)
    ha = _as_array("hartmann", hartmann)
    if np.any(ha <= 0.0):
        raise ValueError("hartmann must be positive")
    return half_width / ha


def cells_per_hartmann_layer(
    cell_size: float | np.ndarray,
    half_width: float,
    hartmann: float | np.ndarray,
) -> float | np.ndarray:
    '''Return the number of uniform wall-normal cells across delta_H.'''

    _require_positive("half_width", half_width)
    dx = _as_array("cell_size", cell_size)
    if np.any(dx <= 0.0):
        raise ValueError("cell_size must be positive")
    return hartmann_layer_thickness(half_width, hartmann) / dx


def classify_hartmann_resolution(
    cell_size: float | np.ndarray,
    half_width: float,
    hartmann: float | np.ndarray,
    target_cells: float = 8.0,
    marginal_fraction: float = 0.5,
) -> str | np.ndarray:
    '''Classify Hartmann-layer resolution as resolved, marginal, or underresolved.'''

    _require_positive("target_cells", target_cells)
    _require_positive("marginal_fraction", marginal_fraction)
    if marginal_fraction >= 1.0:
        raise ValueError("marginal_fraction must be less than 1")
    cpl = cells_per_hartmann_layer(cell_size, half_width, hartmann)
    labels = np.where(
        cpl >= target_cells,
        "resolved",
        np.where(cpl >= target_cells * marginal_fraction, "marginal", "underresolved"),
    )
    if labels.ndim == 0:
        return str(labels.item())
    return labels


def hartmann_mean_over_centerline(hartmann: float | np.ndarray) -> float | np.ndarray:
    '''Return u_mean/u_centerline for insulating Hartmann flow.'''

    ha = _as_array("hartmann", hartmann)
    if np.any(ha < 0.0):
        raise ValueError("hartmann must be non-negative")
    out = np.empty_like(ha, dtype=np.float64)
    small = np.abs(ha) < 1.0e-5
    x = ha[small]
    out[small] = (2.0 / 3.0) + (1.0 / 90.0) * x**2 - (1.0 / 2520.0) * x**4
    xs = ha[~small]
    out[~small] = (1.0 - np.tanh(xs) / xs) / (1.0 - 1.0 / np.cosh(xs))
    return _maybe_scalar(out)


def hartmann_channel_effective_drag(
    half_width: float,
    hartmann: float | np.ndarray,
    nu: float,
) -> float | np.ndarray:
    '''Return alpha_eff matching the insulating Hartmann mean-flow response.

    For a steady channel driven by a uniform streamwise body force per unit mass
    f, the analytic mean velocity is

        u_mean = f a^2 / nu * [1 - tanh(Ha)/Ha] / Ha^2.

    A zero-dimensional reduced balance f = alpha_eff u_mean therefore gives

        alpha_eff = (nu/a^2) Ha^2 / [1 - tanh(Ha)/Ha].

    The Ha -> 0 limit is the Poiseuille mean-flow drag 3 nu/a^2.
    '''

    _require_positive("half_width", half_width)
    _require_positive("nu", nu)
    ha = _as_array("hartmann", hartmann)
    if np.any(ha < 0.0):
        raise ValueError("hartmann must be non-negative")
    factor = _hartmann_mean_force_factor(ha)
    out = (nu / half_width**2) / factor
    return _maybe_scalar(out)


def perpendicular_velocity(velocity: np.ndarray, bhat: np.ndarray) -> np.ndarray:
    '''Return velocity perpendicular to bhat, normalizing bhat if needed.'''

    vel = _as_vector_array("velocity", velocity)
    b = _unit_vector(bhat)
    projection = np.sum(vel * b, axis=-1, keepdims=True)
    return vel - projection * b


def parallel_velocity(velocity: np.ndarray, bhat: np.ndarray) -> np.ndarray:
    '''Return velocity parallel to bhat, normalizing bhat if needed.'''

    vel = _as_vector_array("velocity", velocity)
    b = _unit_vector(bhat)
    projection = np.sum(vel * b, axis=-1, keepdims=True)
    return projection * b


def implicit_drag_update(
    velocity: np.ndarray,
    alpha: float,
    dt: float,
    bhat: np.ndarray,
) -> np.ndarray:
    '''Return backward-Euler update for dv/dt = -alpha v_perp.'''

    _require_nonnegative("alpha", alpha)
    _require_nonnegative("dt", dt)
    vel = _as_vector_array("velocity", velocity)
    vpar = parallel_velocity(vel, bhat)
    vperp = vel - vpar
    return vpar + vperp / (1.0 + alpha * dt)


def drag_power_density(
    velocity: np.ndarray,
    alpha: float,
    bhat: np.ndarray,
    density: float = 1.0,
) -> float | np.ndarray:
    '''Return u dot rho*(-alpha u_perp), which is always non-positive.'''

    _require_nonnegative("alpha", alpha)
    _require_positive("density", density)
    vperp = perpendicular_velocity(velocity, bhat)
    power = -density * alpha * np.sum(vperp * vperp, axis=-1)
    return _maybe_scalar(power)


def _hartmann_mean_force_factor(ha: np.ndarray) -> np.ndarray:
    '''Return [1 - tanh(Ha)/Ha] / Ha^2 with the Ha -> 0 limit.'''

    out = np.empty_like(ha, dtype=np.float64)
    small = np.abs(ha) < 1.0e-5
    x = ha[small]
    out[small] = (1.0 / 3.0) - (2.0 / 15.0) * x**2 + (17.0 / 315.0) * x**4
    xs = ha[~small]
    out[~small] = (1.0 - np.tanh(xs) / xs) / xs**2
    return out


def _as_array(name: str, values: float | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _as_vector_array(name: str, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape == () or arr.shape[-1] != 3:
        raise ValueError(f"{name} must have shape (..., 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _unit_vector(values: np.ndarray) -> np.ndarray:
    vec = np.asarray(values, dtype=np.float64)
    if vec.shape != (3,):
        raise ValueError("bhat must have shape (3,)")
    norm = np.linalg.norm(vec)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("bhat must have positive finite norm")
    return vec / norm


def _require_positive(name: str, value: float) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive")


def _require_nonnegative(name: str, value: float) -> None:
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be non-negative")


def _maybe_scalar(values: np.ndarray) -> float | np.ndarray:
    if values.ndim == 0:
        return float(values)
    return values
