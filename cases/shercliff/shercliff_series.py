"""
Shercliff duct flow — exact series solution for the all-insulating rectangular duct
(Shercliff 1953, Proc. Camb. Phil. Soc. 49, 136).

This is the analytic reference for VALIDATION_CASES.md case S1 (lm-mhd-comsol repo):
square duct, insulating walls (c = 0), Ha = 1323.

Nondimensional problem (duct |x| <= d perpendicular to B, |y| <= 1 along B; axial
velocity u(x,y) driven by unit pressure gradient; b = axial induced field):

    lap(u) + Ha * db/dy = -1        (momentum)
    lap(b) + Ha * du/dy =  0        (induction)

    u = 0 on all four walls (no-slip)
    b = 0 on all four walls (all walls insulating — Shercliff's case)

Method: expand u, b in cos(alpha_k x) modes (alpha_k = (2k+1) pi / (2 d), which
satisfies the insulating side-wall condition term-by-term). The Elsasser variables
p = f + g, m = f - g decouple per mode into constant-coefficient ODEs in y whose
exponents are the Hartmann-layer roots. Exponentials are kept in scaled form
(exp(r*(y -/+ 1)) with |arg| <= 0) so the solution is numerically stable at
arbitrary Ha — no overflow at Ha = 1323 (or 1e6).

Cross-check: solve_fd() is an independent finite-difference referee (same pattern
as cases/hunt/hunt_duct.py but with b = 0 on ALL walls); run with --check to
compare the two at moderate Ha where a uniform FD grid resolves the layers.

Conventions match VALIDATION_CASES.md: profiles vs y/a in [-1, 1]; Hartmann layers
at y = +/-1 (thickness ~ 1/Ha); side layers at x = +/-d (thickness ~ 1/sqrt(Ha)).
Velocity unit is G*a^2/(rho*nu) (unit-pressure-gradient scaling); the contract's
U0 is the reported core velocity u(0,0) in these units.
"""
import argparse
import numpy as np


def _mode_profiles(Ha, alpha, c_k, y):
    """Exact per-mode solution f_k(y) (velocity) and g_k(y) (induced field).

    p = f + g satisfies  p'' + Ha p' - alpha^2 p = -c_k,  p(+/-1) = 0
    m = f - g satisfies  m'' - Ha m' - alpha^2 m = -c_k,  m(+/-1) = 0
    """
    disc = np.sqrt(Ha * Ha + 4.0 * alpha * alpha)
    # roots for p: r1 > 0 small, r2 < 0 ~ -Ha (Hartmann layer at y = -1... y = +1)
    r1 = 0.5 * (-Ha + disc)
    r2 = 0.5 * (-Ha - disc)
    cp = c_k / (alpha * alpha)  # particular solution
    # scaled basis: phi1 = exp(r1*(y-1)) (<=1 since r1>0, y<=1), phi2 = exp(r2*(y+1))
    # BC rows: [phi1(1), phi2(1); phi1(-1), phi2(-1)] * [A, B]^T = [-cp, -cp]^T
    e1m = np.exp(-2.0 * r1)   # phi1(-1)
    e2p = np.exp(2.0 * r2)    # phi2(+1)  (~ exp(-2Ha) -> 0)
    det = 1.0 - e1m * e2p
    A = (-cp * (1.0 - e2p)) / det
    B = (-cp * (1.0 - e1m)) / det
    p = cp + A * np.exp(r1 * (y - 1.0)) + B * np.exp(r2 * (y + 1.0))
    # m(y) = p(-y) by the Ha -> -Ha symmetry of the ODE; use it directly
    m = cp + A * np.exp(r1 * (-y - 1.0)) + B * np.exp(r2 * (-y + 1.0))
    f = 0.5 * (p + m)
    g = 0.5 * (p - m)
    return f, g


def solve_series(Ha, d=1.0, nmodes=4000, x=None, y=None):
    """Exact Shercliff solution on a grid. Returns x, y, u(ny,nx), b(ny,nx)."""
    if x is None:
        x = np.linspace(-d, d, 201)
    if y is None:
        y = np.linspace(-1.0, 1.0, 201)
    u = np.zeros((y.size, x.size))
    b = np.zeros((y.size, x.size))
    for k in range(nmodes):
        alpha = (2 * k + 1) * np.pi / (2.0 * d)
        c_k = 2.0 * (-1.0) ** k / (alpha * d)   # Fourier coeff of RHS "1" on |x|<d
        f, g = _mode_profiles(Ha, alpha, c_k, y)
        ck_x = np.cos(alpha * x)
        u += np.outer(f, ck_x)
        b += np.outer(g, ck_x)
    return x, y, u, b


def core_velocity(Ha, d=1.0, nmodes=4000):
    """u(0,0) — the contract's U0 in units of G*a^2/(rho*nu)."""
    _, _, u, _ = solve_series(Ha, d, nmodes, x=np.array([0.0]), y=np.array([0.0]))
    return float(u[0, 0])


def mean_velocity(Ha, d=1.0, nmodes=4000, n=801):
    """Cross-section mean of u (for flow-rate normalizations)."""
    x = np.linspace(-d, d, n)
    y = np.linspace(-1.0, 1.0, n)
    _, _, u, _ = solve_series(Ha, d, nmodes, x=x, y=y)
    return float(np.trapezoid(np.trapezoid(u, x, axis=1), y) / (2.0 * d * 2.0))


def solve_fd(Ha, d=1.0, nx=201, ny=201):
    """Independent FD referee: coupled equations, b = 0 on all walls."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    x = np.linspace(-d, d, nx)
    y = np.linspace(-1.0, 1.0, ny)
    dx, dy = x[1] - x[0], y[1] - y[0]
    N = nx * ny
    idx = lambda i, j: j * nx + i
    A = sp.lil_matrix((2 * N, 2 * N))
    rhs = np.zeros(2 * N)
    for j in range(ny):
        for i in range(nx):
            uId, bId = idx(i, j), N + idx(i, j)
            wall = i in (0, nx - 1) or j in (0, ny - 1)
            if wall:
                A[uId, uId] = 1.0
                A[bId, bId] = 1.0
                continue
            for (rId, oId) in ((uId, N), (bId, -N)):
                A[rId, rId] = -2.0 / dx**2 - 2.0 / dy**2
                A[rId, rId - 1] = 1.0 / dx**2
                A[rId, rId + 1] = 1.0 / dx**2
                A[rId, rId - nx] = 1.0 / dy**2
                A[rId, rId + nx] = 1.0 / dy**2
                A[rId, rId + oId + nx] = Ha / (2.0 * dy)
                A[rId, rId + oId - nx] = -Ha / (2.0 * dy)
            rhs[uId] = -1.0
    z = spla.spsolve(A.tocsr(), rhs)
    return x, y, z[:N].reshape(ny, nx), z[N:].reshape(ny, nx)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--Ha", type=float, default=1323.0)
    ap.add_argument("--d", type=float, default=1.0, help="aspect: side-wall half-width / a")
    ap.add_argument("--nmodes", type=int, default=4000)
    ap.add_argument("--check", action="store_true", help="FD cross-check at this Ha")
    ap.add_argument("--profile", metavar="CSV", help="write centerline profiles to CSV")
    args = ap.parse_args()

    U0 = core_velocity(args.Ha, args.d, args.nmodes)
    Um = mean_velocity(args.Ha, args.d, args.nmodes)
    print(f"Ha = {args.Ha:g}, square duct d = {args.d:g}, {args.nmodes} modes")
    print(f"  U0  (core velocity u(0,0), units G*a^2/(rho*nu)) = {U0:.8e}")
    print(f"  Um  (mean velocity, same units)                  = {Um:.8e}")
    print(f"  U0 * Ha = {U0 * args.Ha:.6f}   (Hartmann-flow limit -> 1)")

    if args.check:
        nx = ny = 241
        xs, ys, u_fd, b_fd = solve_fd(args.Ha, args.d, nx, ny)
        _, _, u_se, b_se = solve_series(args.Ha, args.d, args.nmodes, x=xs, y=ys)
        err = np.max(np.abs(u_fd - u_se)) / np.max(np.abs(u_se))
        print(f"  FD referee ({nx}x{ny}): max |u_fd - u_series| / max|u| = {err:.3e}")

    if args.profile:
        ys = np.linspace(-1.0, 1.0, 1001)
        xs = np.linspace(-args.d, args.d, 1001)
        _, _, u_y, _ = solve_series(args.Ha, args.d, args.nmodes,
                                    x=np.array([0.0]), y=ys)
        _, _, u_x, _ = solve_series(args.Ha, args.d, args.nmodes,
                                    x=xs, y=np.array([0.0]))
        with open(args.profile, "w") as fh:
            fh.write(f"# Shercliff series, Ha={args.Ha:g}, d={args.d:g}, "
                     f"nmodes={args.nmodes}, U0={U0:.8e}\n")
            fh.write("# col1: coord, col2: u(0,y) [Hartmann cut], col3: u(x,0) [side cut]\n")
            for i in range(ys.size):
                fh.write(f"{ys[i]:.6f} {u_y[i,0]:.8e} {u_x[0,i]:.8e}\n")
        print(f"  wrote {args.profile}")
