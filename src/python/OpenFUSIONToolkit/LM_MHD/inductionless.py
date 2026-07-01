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
import scipy.sparse as sp
import scipy.sparse.linalg as spla


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


def structured_gradient_2d(
    scalar: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> np.ndarray:
    '''Return grad(scalar) on a structured x-z grid as a 3-vector field.'''

    values, xcoord, zcoord = _structured_scalar("scalar", scalar, x, z)
    comp_x, comp_z = _check_components(components)
    gradient = np.zeros(values.shape + (3,), dtype=np.float64)
    gradient[..., comp_x] = np.gradient(values, xcoord, axis=1, edge_order=2)
    gradient[..., comp_z] = np.gradient(values, zcoord, axis=0, edge_order=2)
    return gradient


def motional_electric_field(velocity: np.ndarray, magnetic_field: np.ndarray) -> np.ndarray:
    '''Return the motional term u x B used in inductionless Ohm's law.'''

    vel = _as_vector_array("velocity", velocity)
    bfield = _broadcast_vector("magnetic_field", magnetic_field, vel.shape)
    return np.cross(vel, bfield)


def potential_source_from_motional_emf(
    motional_emf: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> np.ndarray:
    '''Return source div(u x B) for constant-conductivity potential solves.'''

    emf = _structured_vector("motional_emf", motional_emf, x, z)
    comp_x, comp_z = _check_components(components)
    xcoord, zcoord = _structured_coordinates(x, z)
    d_ex_dx = np.gradient(emf[..., comp_x], xcoord, axis=1, edge_order=2)
    d_ez_dz = np.gradient(emf[..., comp_z], zcoord, axis=0, edge_order=2)
    return d_ex_dx + d_ez_dz


def insulating_wall_normal_gradient(
    motional_emf: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> dict[str, np.ndarray]:
    '''Return d(phi)/dn = (u x B).n for insulating walls.

    The returned arrays are outward-normal derivatives for the rectangle sides
    expected by `solve_potential_neumann_2d`.
    '''

    emf = _structured_vector("motional_emf", motional_emf, x, z)
    comp_x, comp_z = _check_components(components)
    return {
        "x_min": -emf[:, 0, comp_x].copy(),
        "x_max": emf[:, -1, comp_x].copy(),
        "z_min": -emf[0, :, comp_z].copy(),
        "z_max": emf[-1, :, comp_z].copy(),
    }


def reconstruct_inductionless_current_2d(
    potential: np.ndarray,
    velocity: np.ndarray,
    magnetic_field: np.ndarray,
    sigma: float | np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> np.ndarray:
    '''Return J = sigma(-grad(phi) + u x B) on a structured x-z grid.'''

    grad_phi = structured_gradient_2d(potential, x, z, components=components)
    return ohms_law_current_density(sigma, velocity, magnetic_field, grad_phi)


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


def solve_potential_dirichlet_2d(
    source: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    boundary_values: float | np.ndarray = 0.0,
) -> np.ndarray:
    '''Solve laplacian(phi) = source on a uniform rectangular grid.

    This is a finite-difference reference solve for manufactured tests of the
    future inductionless electric-potential equation. It uses Dirichlet boundary
    values on all four sides. Insulating Neumann and conducting-wall Robin
    closures are deliberately not represented by this helper.
    '''

    src, xcoord, zcoord = _structured_scalar("source", source, x, z)
    hx = _uniform_spacing("x", xcoord)
    hz = _uniform_spacing("z", zcoord)
    boundary = _boundary_values(boundary_values, src.shape)
    nz, nx = src.shape
    nxi = nx - 2
    nzi = nz - 2
    if nxi < 1 or nzi < 1:
        raise ValueError("source grid must include at least one interior point")

    rows = []
    cols = []
    data = []
    rhs = np.empty(nxi * nzi, dtype=np.float64)
    inv_hx2 = 1.0 / hx**2
    inv_hz2 = 1.0 / hz**2

    def row_index(iz: int, ix: int) -> int:
        return iz * nxi + ix

    for iz in range(nzi):
        for ix in range(nxi):
            row = row_index(iz, ix)
            rows.append(row)
            cols.append(row)
            data.append(-2.0 * (inv_hx2 + inv_hz2))
            value = src[iz + 1, ix + 1]

            if ix > 0:
                rows.append(row)
                cols.append(row_index(iz, ix - 1))
                data.append(inv_hx2)
            else:
                value -= boundary[iz + 1, 0] * inv_hx2
            if ix < nxi - 1:
                rows.append(row)
                cols.append(row_index(iz, ix + 1))
                data.append(inv_hx2)
            else:
                value -= boundary[iz + 1, -1] * inv_hx2

            if iz > 0:
                rows.append(row)
                cols.append(row_index(iz - 1, ix))
                data.append(inv_hz2)
            else:
                value -= boundary[0, ix + 1] * inv_hz2
            if iz < nzi - 1:
                rows.append(row)
                cols.append(row_index(iz + 1, ix))
                data.append(inv_hz2)
            else:
                value -= boundary[-1, ix + 1] * inv_hz2
            rhs[row] = value

    matrix = sp.csr_matrix((data, (rows, cols)), shape=(nxi * nzi, nxi * nzi))
    interior = spla.spsolve(matrix, rhs)
    phi = boundary.copy()
    phi[1:-1, 1:-1] = interior.reshape(nzi, nxi)
    return phi


def solve_variable_conductivity_dirichlet_2d(
    conductivity: np.ndarray,
    source: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    boundary_values: float | np.ndarray = 0.0,
) -> np.ndarray:
    '''Solve div(sigma grad(phi)) = source with Dirichlet walls.

    This uniform-grid reference uses harmonic face averages for sigma. It is
    meant for manufactured tests of regional/material-coefficient handling, not
    for production blanket solves.
    '''

    sigma, xcoord, zcoord = _structured_scalar("conductivity", conductivity, x, z)
    if np.any(sigma <= 0.0):
        raise ValueError("conductivity must be positive")
    src, _, _ = _structured_scalar("source", source, x, z)
    hx = _uniform_spacing("x", xcoord)
    hz = _uniform_spacing("z", zcoord)
    boundary = _boundary_values(boundary_values, src.shape)
    nz, nx = src.shape
    nxi = nx - 2
    nzi = nz - 2
    if nxi < 1 or nzi < 1:
        raise ValueError("source grid must include at least one interior point")

    rows = []
    cols = []
    data = []
    rhs = np.empty(nxi * nzi, dtype=np.float64)
    inv_hx2 = 1.0 / hx**2
    inv_hz2 = 1.0 / hz**2

    def row_index(iz: int, ix: int) -> int:
        return iz * nxi + ix

    for iz in range(nzi):
        for ix in range(nxi):
            grid_z = iz + 1
            grid_x = ix + 1
            row = row_index(iz, ix)
            sigma_center = sigma[grid_z, grid_x]
            sigma_e = _harmonic_mean(sigma_center, sigma[grid_z, grid_x + 1])
            sigma_w = _harmonic_mean(sigma_center, sigma[grid_z, grid_x - 1])
            sigma_n = _harmonic_mean(sigma_center, sigma[grid_z + 1, grid_x])
            sigma_s = _harmonic_mean(sigma_center, sigma[grid_z - 1, grid_x])

            rows.append(row)
            cols.append(row)
            data.append(-(sigma_e + sigma_w) * inv_hx2 - (sigma_n + sigma_s) * inv_hz2)
            value = src[grid_z, grid_x]

            west_coeff = sigma_w * inv_hx2
            east_coeff = sigma_e * inv_hx2
            south_coeff = sigma_s * inv_hz2
            north_coeff = sigma_n * inv_hz2

            if ix > 0:
                rows.append(row)
                cols.append(row_index(iz, ix - 1))
                data.append(west_coeff)
            else:
                value -= west_coeff * boundary[grid_z, 0]
            if ix < nxi - 1:
                rows.append(row)
                cols.append(row_index(iz, ix + 1))
                data.append(east_coeff)
            else:
                value -= east_coeff * boundary[grid_z, -1]

            if iz > 0:
                rows.append(row)
                cols.append(row_index(iz - 1, ix))
                data.append(south_coeff)
            else:
                value -= south_coeff * boundary[0, grid_x]
            if iz < nzi - 1:
                rows.append(row)
                cols.append(row_index(iz + 1, ix))
                data.append(north_coeff)
            else:
                value -= north_coeff * boundary[-1, grid_x]
            rhs[row] = value

    matrix = sp.csr_matrix((data, (rows, cols)), shape=(nxi * nzi, nxi * nzi))
    interior = spla.spsolve(matrix, rhs)
    phi = boundary.copy()
    phi[1:-1, 1:-1] = interior.reshape(nzi, nxi)
    return phi


def solve_potential_neumann_2d(
    source: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    mean_value: float = 0.0,
    normal_gradient: dict[str, np.ndarray | float] | None = None,
    compatibility_tolerance: float = 1.0e-9,
) -> np.ndarray:
    '''Solve laplacian(phi) = source with Neumann wall data.

    The solve uses a reflected-ghost finite-difference stencil on a uniform
    rectangular grid and fixes the nullspace by enforcing the requested mean
    value. `normal_gradient` gives outward-normal d(phi)/dn data on the four
    sides: `x_min`, `x_max`, `z_min`, and `z_max`. Missing sides default to
    zero. Robin conducting-wall closure is not represented by this helper.
    '''

    src, xcoord, zcoord = _structured_scalar("source", source, x, z)
    _require_finite("mean_value", mean_value)
    _require_positive("compatibility_tolerance", compatibility_tolerance)
    hx = _uniform_spacing("x", xcoord)
    hz = _uniform_spacing("z", zcoord)
    nz, nx = src.shape
    gradients = _normal_gradient_values(normal_gradient, nx, nz)
    rhs_grid = src.copy()
    size = nx * nz
    rows = []
    cols = []
    data = []
    inv_hx2 = 1.0 / hx**2
    inv_hz2 = 1.0 / hz**2

    def row_index(iz: int, ix: int) -> int:
        return iz * nx + ix

    for iz in range(nz):
        for ix in range(nx):
            row = row_index(iz, ix)
            if ix == 0:
                rows.extend((row, row))
                cols.extend((row_index(iz, 0), row_index(iz, 1)))
                data.extend((-2.0 * inv_hx2, 2.0 * inv_hx2))
                rhs_grid[iz, ix] -= 2.0 * gradients["x_min"][iz] / hx
            elif ix == nx - 1:
                rows.extend((row, row))
                cols.extend((row_index(iz, nx - 1), row_index(iz, nx - 2)))
                data.extend((-2.0 * inv_hx2, 2.0 * inv_hx2))
                rhs_grid[iz, ix] -= 2.0 * gradients["x_max"][iz] / hx
            else:
                rows.extend((row, row, row))
                cols.extend((row, row_index(iz, ix - 1), row_index(iz, ix + 1)))
                data.extend((-2.0 * inv_hx2, inv_hx2, inv_hx2))

            if iz == 0:
                rows.extend((row, row))
                cols.extend((row_index(0, ix), row_index(1, ix)))
                data.extend((-2.0 * inv_hz2, 2.0 * inv_hz2))
                rhs_grid[iz, ix] -= 2.0 * gradients["z_min"][ix] / hz
            elif iz == nz - 1:
                rows.extend((row, row))
                cols.extend((row_index(nz - 1, ix), row_index(nz - 2, ix)))
                data.extend((-2.0 * inv_hz2, 2.0 * inv_hz2))
                rhs_grid[iz, ix] -= 2.0 * gradients["z_max"][ix] / hz
            else:
                rows.extend((row, row, row))
                cols.extend((row, row_index(iz - 1, ix), row_index(iz + 1, ix)))
                data.extend((-2.0 * inv_hz2, inv_hz2, inv_hz2))

    matrix = sp.csr_matrix((data, (rows, cols)), shape=(size, size))
    constraint = sp.csr_matrix(np.ones((1, size), dtype=np.float64))
    augmented = sp.bmat([[matrix, constraint.T], [constraint, None]], format="csr")
    rhs = np.concatenate((rhs_grid.ravel(), np.array([float(mean_value) * size])))
    solution = spla.spsolve(augmented, rhs)
    gauge_multiplier = float(solution[-1])
    rhs_scale = max(1.0, float(np.max(np.abs(rhs_grid))))
    if abs(gauge_multiplier) > compatibility_tolerance * rhs_scale:
        raise ValueError("Neumann source and normal_gradient data appear incompatible")
    return solution[:size].reshape(nz, nx)


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


def _structured_vector(
    name: str,
    values: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    arr = _as_vector_array(name, values)
    if arr.ndim != 3:
        raise ValueError(f"{name} must have shape (nz, nx, 3)")
    xcoord, zcoord = _structured_coordinates(x, z)
    if arr.shape[:2] != (zcoord.size, xcoord.size):
        raise ValueError(f"{name} shape must match z and x coordinates")
    return arr


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


def _boundary_values(values: float | np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    if np.isscalar(values):
        out = np.full(shape, float(values), dtype=np.float64)
    else:
        out = np.asarray(values, dtype=np.float64)
        if out.shape != shape:
            raise ValueError("boundary_values array must match source shape")
        out = out.copy()
    if not np.all(np.isfinite(out)):
        raise ValueError("boundary_values must contain only finite values")
    return out


def _uniform_spacing(name: str, coord: np.ndarray) -> float:
    diffs = np.diff(coord)
    spacing = float(diffs[0])
    if not np.allclose(diffs, spacing, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError(f"{name} must be uniformly spaced for this reference solve")
    return spacing


def _harmonic_mean(left: float, right: float) -> float:
    return float(2.0 * left * right / (left + right))


def _normal_gradient_values(
    values: dict[str, np.ndarray | float] | None,
    nx: int,
    nz: int,
) -> dict[str, np.ndarray]:
    gradients = {
        "x_min": np.zeros(nz, dtype=np.float64),
        "x_max": np.zeros(nz, dtype=np.float64),
        "z_min": np.zeros(nx, dtype=np.float64),
        "z_max": np.zeros(nx, dtype=np.float64),
    }
    if values is None:
        return gradients
    allowed = set(gradients)
    extra = set(values) - allowed
    if extra:
        raise ValueError(f"unknown normal_gradient side(s): {sorted(extra)}")
    for side, target in gradients.items():
        if side not in values:
            continue
        raw = np.asarray(values[side], dtype=np.float64)
        if raw.shape == ():
            out = np.full(target.shape, float(raw), dtype=np.float64)
        else:
            if raw.shape != target.shape:
                raise ValueError(f"normal_gradient['{side}'] must have shape {target.shape}")
            out = raw.copy()
        if not np.all(np.isfinite(out)):
            raise ValueError(f"normal_gradient['{side}'] must contain only finite values")
        gradients[side] = out
    return gradients


def _require_positive(name: str, value: float) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive")


def _require_finite(name: str, value: float) -> None:
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")


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
