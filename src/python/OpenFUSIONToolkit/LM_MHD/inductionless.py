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


def weighted_potential_source_from_motional_emf(
    conductivity: np.ndarray,
    motional_emf: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> np.ndarray:
    '''Return div(sigma * (u x B)) for variable-conductivity solves.'''

    sigma, xcoord, zcoord = _structured_scalar("conductivity", conductivity, x, z)
    if np.any(sigma <= 0.0):
        raise ValueError("conductivity must be positive")
    emf = _structured_vector("motional_emf", motional_emf, x, z)
    comp_x, comp_z = _check_components(components)
    weighted_x = sigma * emf[..., comp_x]
    weighted_z = sigma * emf[..., comp_z]
    d_ex_dx = np.gradient(weighted_x, xcoord, axis=1, edge_order=2)
    d_ez_dz = np.gradient(weighted_z, zcoord, axis=0, edge_order=2)
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


def insulating_wall_normal_flux(
    conductivity: np.ndarray,
    motional_emf: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    components: tuple[int, int] = (0, 2),
) -> dict[str, np.ndarray]:
    '''Return sigma*(u x B).n wall flux for insulating current closure.'''

    sigma, _, _ = _structured_scalar("conductivity", conductivity, x, z)
    if np.any(sigma <= 0.0):
        raise ValueError("conductivity must be positive")
    emf = _structured_vector("motional_emf", motional_emf, x, z)
    comp_x, comp_z = _check_components(components)
    return {
        "x_min": -(sigma[:, 0] * emf[:, 0, comp_x]).copy(),
        "x_max": (sigma[:, -1] * emf[:, -1, comp_x]).copy(),
        "z_min": -(sigma[0, :] * emf[0, :, comp_z]).copy(),
        "z_max": (sigma[-1, :] * emf[-1, :, comp_z]).copy(),
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


def axisymmetric_electric_field_from_potential(
    potential: np.ndarray,
    r: np.ndarray,
    z: np.ndarray,
    toroidal_electric_field: float | np.ndarray = 0.0,
) -> np.ndarray:
    '''Return axisymmetric electric field E = -grad(phi) + E_phi e_phi.

    Components are ordered `(R, phi, Z)`.
    '''

    potential_values, rcoord, zcoord = _structured_scalar("potential", potential, r, z)
    electric_field = np.zeros(potential_values.shape + (3,), dtype=np.float64)
    electric_field[..., 0] = -np.gradient(potential_values, rcoord, axis=1, edge_order=2)
    electric_field[..., 1] = _broadcast_scalar_field(
        "toroidal_electric_field",
        toroidal_electric_field,
        potential_values.shape,
    )
    electric_field[..., 2] = -np.gradient(potential_values, zcoord, axis=0, edge_order=2)
    return electric_field


def reconstruct_axisymmetric_inductionless_current(
    potential: np.ndarray,
    velocity: np.ndarray,
    magnetic_field: np.ndarray,
    sigma: float | np.ndarray,
    r: np.ndarray,
    z: np.ndarray,
    toroidal_electric_field: float | np.ndarray = 0.0,
) -> np.ndarray:
    '''Return axisymmetric J = sigma(-grad(phi) + u x B + E_phi e_phi).

    Components are ordered `(R, phi, Z)`. The scalar potential has only
    poloidal gradients under axisymmetry; `toroidal_electric_field` represents
    an externally imposed loop-voltage/toroidal-electric-field contribution.
    '''

    electric_field = axisymmetric_electric_field_from_potential(
        potential,
        r,
        z,
        toroidal_electric_field=toroidal_electric_field,
    )
    vel = _broadcast_vector("velocity", velocity, electric_field.shape)
    bfield = _broadcast_vector("magnetic_field", magnetic_field, electric_field.shape)
    sigma_arr = _positive_broadcast("sigma", sigma, electric_field.shape[:-1])
    return sigma_arr[..., None] * (electric_field + np.cross(vel, bfield))


def toroidal_electric_field_from_loop_voltage(
    loop_voltage: float | np.ndarray,
    r: float | np.ndarray,
) -> float | np.ndarray:
    '''Return axisymmetric E_phi = V_loop / (2*pi*R).'''

    radius = np.asarray(r, dtype=np.float64)
    if not np.all(np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("r must contain positive finite values")
    voltage = np.asarray(loop_voltage, dtype=np.float64)
    if not np.all(np.isfinite(voltage)):
        raise ValueError("loop_voltage must contain only finite values")
    field = voltage / (2.0 * np.pi * radius)
    return _maybe_scalar(field)


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


def axisymmetric_current_from_flux_function(
    flux_function: np.ndarray,
    r: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    '''Return axisymmetric divergence-free poloidal current from a flux function.

    Components are `(R, phi, Z)` with
    `J_R = (1/R) dpsi/dZ`, `J_Z = -(1/R) dpsi/dR`, and `J_phi = 0`.
    This makes `(1/R) d(R J_R)/dR + dJ_Z/dZ = 0` up to the structured-grid
    differentiation error.
    '''

    psi, rcoord, zcoord = _structured_scalar("flux_function", flux_function, r, z)
    if np.any(rcoord <= 0.0):
        raise ValueError("r must be positive for axisymmetric current diagnostics")
    r_grid = rcoord[None, :]
    current = np.zeros(psi.shape + (3,), dtype=np.float64)
    current[..., 0] = np.gradient(psi, zcoord, axis=0, edge_order=2) / r_grid
    current[..., 2] = -np.gradient(psi, rcoord, axis=1, edge_order=2) / r_grid
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


def axisymmetric_charge_conservation_residual(
    current_density: np.ndarray,
    r: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    '''Return axisymmetric residual div(J) = (1/R)d(RJ_R)/dR + dJ_Z/dZ.'''

    current = _as_vector_array("current_density", current_density)
    if current.ndim != 3:
        raise ValueError("current_density must have shape (nz, nr, 3)")
    rcoord, zcoord = _structured_coordinates(r, z)
    if np.any(rcoord <= 0.0):
        raise ValueError("r must be positive for axisymmetric divergence")
    if current.shape[:2] != (zcoord.size, rcoord.size):
        raise ValueError("current_density shape must match z and r coordinates")
    r_grid = rcoord[None, :]
    d_rjr_dr = np.gradient(r_grid * current[..., 0], rcoord, axis=1, edge_order=2)
    d_jz_dz = np.gradient(current[..., 2], zcoord, axis=0, edge_order=2)
    return d_rjr_dr / r_grid + d_jz_dz


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


def axisymmetric_wall_normal_current_extrema(current_density: np.ndarray) -> dict[str, float]:
    '''Return max absolute normal current on rectangular R-Z grid boundaries.'''

    current = _as_vector_array("current_density", current_density)
    if current.ndim != 3:
        raise ValueError("current_density must have shape (nz, nr, 3)")
    return {
        "r_min": float(np.max(np.abs(current[:, 0, 0]))),
        "r_max": float(np.max(np.abs(current[:, -1, 0]))),
        "z_min": float(np.max(np.abs(current[0, :, 2]))),
        "z_max": float(np.max(np.abs(current[-1, :, 2]))),
    }


def axisymmetric_boundary_current_integrals(current_density: np.ndarray, r: np.ndarray, z: np.ndarray) -> dict[str, float]:
    '''Return signed outward boundary-current integrals on an R-Z grid.

    Components are `(R, phi, Z)`. The returned values have units of amperes.
    Radial surfaces use `dS = 2*pi*R*dZ`; horizontal surfaces use
    `dS = 2*pi*R*dR`.
    '''

    current = _as_vector_array("current_density", current_density)
    if current.ndim != 3:
        raise ValueError("current_density must have shape (nz, nr, 3)")
    rcoord, zcoord = _structured_coordinates(r, z)
    if np.any(rcoord <= 0.0):
        raise ValueError("r must be positive for axisymmetric boundary integrals")
    if current.shape[:2] != (zcoord.size, rcoord.size):
        raise ValueError("current_density shape must match z and r coordinates")
    dr = _node_widths(rcoord)
    dz = _node_widths(zcoord)
    return {
        "r_min": float(np.sum(-current[:, 0, 0] * (2.0 * np.pi * rcoord[0] * dz))),
        "r_max": float(np.sum(current[:, -1, 0] * (2.0 * np.pi * rcoord[-1] * dz))),
        "z_min": float(np.sum(-current[0, :, 2] * (2.0 * np.pi * rcoord * dr))),
        "z_max": float(np.sum(current[-1, :, 2] * (2.0 * np.pi * rcoord * dr))),
    }


def axisymmetric_net_boundary_current(current_density: np.ndarray, r: np.ndarray, z: np.ndarray) -> float:
    '''Return net signed outward current through the R-Z boundary.'''

    boundary_current = axisymmetric_boundary_current_integrals(current_density, r, z)
    return float(sum(boundary_current.values()))


def axisymmetric_cell_volumes(r: np.ndarray, z: np.ndarray) -> np.ndarray:
    '''Return midpoint-control-volume weights dV = 2*pi*R*dR*dZ on an R-Z grid.'''

    rcoord, zcoord = _structured_coordinates(r, z)
    if np.any(rcoord <= 0.0):
        raise ValueError("r must be positive for axisymmetric volumes")
    dr = _node_widths(rcoord)
    dz = _node_widths(zcoord)
    return 2.0 * np.pi * rcoord[None, :] * dz[:, None] * dr[None, :]


def mechanical_power_density(force_density: np.ndarray, velocity: np.ndarray) -> float | np.ndarray:
    '''Return local mechanical power density p = f.u.'''

    force = _as_vector_array("force_density", force_density)
    vel = _broadcast_vector("velocity", velocity, force.shape)
    power = np.sum(force * vel, axis=-1)
    return _maybe_scalar(power)


def electric_power_density(current_density: np.ndarray, electric_field: np.ndarray) -> float | np.ndarray:
    '''Return local electrical power density p = J.E.'''

    current = _as_vector_array("current_density", current_density)
    electric = _broadcast_vector("electric_field", electric_field, current.shape)
    power = np.sum(current * electric, axis=-1)
    return _maybe_scalar(power)


def inductionless_power_balance_residual(
    current_density: np.ndarray,
    velocity: np.ndarray,
    magnetic_field: np.ndarray,
    sigma: float | np.ndarray,
    electric_field: np.ndarray,
) -> float | np.ndarray:
    '''Return q_J + (J x B).u - J.E for inductionless Ohm's law.

    For `J = sigma(E + u x B)`, `q_J = |J|^2/sigma`, and
    `f = J x B`, this residual should vanish pointwise up to roundoff.
    '''

    current = _as_vector_array("current_density", current_density)
    force = lorentz_force_density(current, magnetic_field)
    joule = joule_heating_density(current, sigma)
    mechanical = mechanical_power_density(force, velocity)
    electric = electric_power_density(current, electric_field)
    residual = joule + mechanical - electric
    return _maybe_scalar(np.asarray(residual))


def axisymmetric_torque_density(
    force_density: np.ndarray,
    r: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    '''Return torque-about-axis density tau = R*f_phi on an R-Z grid.'''

    force = _as_vector_array("force_density", force_density)
    if force.ndim != 3:
        raise ValueError("force_density must have shape (nz, nr, 3)")
    rcoord, zcoord = _structured_coordinates(r, z)
    if np.any(rcoord <= 0.0):
        raise ValueError("r must be positive for axisymmetric torque diagnostics")
    if force.shape[:2] != (zcoord.size, rcoord.size):
        raise ValueError("force_density shape must match z and r coordinates")
    return rcoord[None, :] * force[..., 1]


def axisymmetric_volume_integral(density: np.ndarray, r: np.ndarray, z: np.ndarray) -> float:
    '''Return integral density dV over an axisymmetric R-Z grid.'''

    values, rcoord, zcoord = _structured_scalar("density", density, r, z)
    return float(np.sum(values * axisymmetric_cell_volumes(rcoord, zcoord)))


def axisymmetric_integrated_power(power_density: np.ndarray, r: np.ndarray, z: np.ndarray) -> float:
    '''Return volume-integrated power from a scalar power density.'''

    return axisymmetric_volume_integral(power_density, r, z)


def axisymmetric_integrated_torque(force_density: np.ndarray, r: np.ndarray, z: np.ndarray) -> float:
    '''Return torque about the symmetry axis from force density.'''

    return axisymmetric_volume_integral(axisymmetric_torque_density(force_density, r, z), r, z)


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


def solve_variable_conductivity_neumann_2d(
    conductivity: np.ndarray,
    source: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    mean_value: float = 0.0,
    normal_flux: dict[str, np.ndarray | float] | None = None,
    compatibility_tolerance: float = 1.0e-9,
) -> np.ndarray:
    '''Solve div(sigma grad(phi)) = source with Neumann flux data.

    `normal_flux` gives outward sigma*d(phi)/dn on `x_min`, `x_max`,
    `z_min`, and `z_max`. Missing sides default to zero. Harmonic face
    averages are used for interior conductivity.
    '''

    sigma, xcoord, zcoord = _structured_scalar("conductivity", conductivity, x, z)
    if np.any(sigma <= 0.0):
        raise ValueError("conductivity must be positive")
    src, _, _ = _structured_scalar("source", source, x, z)
    _require_finite("mean_value", mean_value)
    _require_positive("compatibility_tolerance", compatibility_tolerance)
    hx = _uniform_spacing("x", xcoord)
    hz = _uniform_spacing("z", zcoord)
    nz, nx = src.shape
    fluxes = _normal_flux_values(normal_flux, nx, nz)
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
            sigma_center = sigma[iz, ix]

            if ix == 0:
                sigma_e = _harmonic_mean(sigma_center, sigma[iz, ix + 1])
                rows.extend((row, row))
                cols.extend((row_index(iz, ix), row_index(iz, ix + 1)))
                data.extend((-2.0 * sigma_e * inv_hx2, 2.0 * sigma_e * inv_hx2))
                rhs_grid[iz, ix] -= 2.0 * fluxes["x_min"][iz] / hx
            elif ix == nx - 1:
                sigma_w = _harmonic_mean(sigma_center, sigma[iz, ix - 1])
                rows.extend((row, row))
                cols.extend((row_index(iz, ix), row_index(iz, ix - 1)))
                data.extend((-2.0 * sigma_w * inv_hx2, 2.0 * sigma_w * inv_hx2))
                rhs_grid[iz, ix] -= 2.0 * fluxes["x_max"][iz] / hx
            else:
                sigma_e = _harmonic_mean(sigma_center, sigma[iz, ix + 1])
                sigma_w = _harmonic_mean(sigma_center, sigma[iz, ix - 1])
                rows.extend((row, row, row))
                cols.extend((row, row_index(iz, ix - 1), row_index(iz, ix + 1)))
                data.extend((-(sigma_e + sigma_w) * inv_hx2, sigma_w * inv_hx2, sigma_e * inv_hx2))

            if iz == 0:
                sigma_n = _harmonic_mean(sigma_center, sigma[iz + 1, ix])
                rows.extend((row, row))
                cols.extend((row_index(iz, ix), row_index(iz + 1, ix)))
                data.extend((-2.0 * sigma_n * inv_hz2, 2.0 * sigma_n * inv_hz2))
                rhs_grid[iz, ix] -= 2.0 * fluxes["z_min"][ix] / hz
            elif iz == nz - 1:
                sigma_s = _harmonic_mean(sigma_center, sigma[iz - 1, ix])
                rows.extend((row, row))
                cols.extend((row_index(iz, ix), row_index(iz - 1, ix)))
                data.extend((-2.0 * sigma_s * inv_hz2, 2.0 * sigma_s * inv_hz2))
                rhs_grid[iz, ix] -= 2.0 * fluxes["z_max"][ix] / hz
            else:
                sigma_n = _harmonic_mean(sigma_center, sigma[iz + 1, ix])
                sigma_s = _harmonic_mean(sigma_center, sigma[iz - 1, ix])
                rows.extend((row, row, row))
                cols.extend((row, row_index(iz - 1, ix), row_index(iz + 1, ix)))
                data.extend((-(sigma_n + sigma_s) * inv_hz2, sigma_s * inv_hz2, sigma_n * inv_hz2))

    matrix = sp.csr_matrix((data, (rows, cols)), shape=(size, size))
    constraint = sp.csr_matrix(np.ones((1, size), dtype=np.float64))
    augmented = sp.bmat([[matrix, constraint.T], [constraint, None]], format="csr")
    rhs = np.concatenate((rhs_grid.ravel(), np.array([float(mean_value) * size])))
    solution = spla.spsolve(augmented, rhs)
    gauge_multiplier = float(solution[-1])
    rhs_scale = max(1.0, float(np.max(np.abs(rhs_grid))))
    if abs(gauge_multiplier) > compatibility_tolerance * rhs_scale:
        raise ValueError("Neumann source and normal_flux data appear incompatible")
    return solution[:size].reshape(nz, nx)


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


def _broadcast_scalar_field(name: str, values: float | np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    try:
        arr = np.broadcast_to(arr, target_shape)
    except ValueError as exc:
        raise ValueError(f"{name} must broadcast to {target_shape}") from exc
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
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


def _node_widths(coord: np.ndarray) -> np.ndarray:
    widths = np.empty_like(coord, dtype=np.float64)
    widths[1:-1] = 0.5 * (coord[2:] - coord[:-2])
    widths[0] = 0.5 * (coord[1] - coord[0])
    widths[-1] = 0.5 * (coord[-1] - coord[-2])
    return widths


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


def _normal_flux_values(
    values: dict[str, np.ndarray | float] | None,
    nx: int,
    nz: int,
) -> dict[str, np.ndarray]:
    fluxes = {
        "x_min": np.zeros(nz, dtype=np.float64),
        "x_max": np.zeros(nz, dtype=np.float64),
        "z_min": np.zeros(nx, dtype=np.float64),
        "z_max": np.zeros(nx, dtype=np.float64),
    }
    if values is None:
        return fluxes
    allowed = set(fluxes)
    extra = set(values) - allowed
    if extra:
        raise ValueError(f"unknown normal_flux side(s): {sorted(extra)}")
    for side, target in fluxes.items():
        if side not in values:
            continue
        raw = np.asarray(values[side], dtype=np.float64)
        if raw.shape == ():
            out = np.full(target.shape, float(raw), dtype=np.float64)
        else:
            if raw.shape != target.shape:
                raise ValueError(f"normal_flux['{side}'] must have shape {target.shape}")
            out = raw.copy()
        if not np.all(np.isfinite(out)):
            raise ValueError(f"normal_flux['{side}'] must contain only finite values")
        fluxes[side] = out
    return fluxes


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
