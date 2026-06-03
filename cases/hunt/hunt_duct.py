"""
Hunt's duct flow — finite-difference reference for the conducting-Hartmann-wall,
insulating-side-wall rectangular duct (Hunt 1965, J. Fluid Mech. 21, 577).

This is a robust *semi-analytic* reference: it solves the standard inductionless
(low-Rm) coupled MHD duct equations directly on a grid, so it unambiguously
exhibits the SIGNATURE physics — high-velocity side jets near the insulating side
walls and a near-uniform core — which distinguishes Hunt's flow from the insulating
Hartmann/Shercliff case. (The closed-form series is attempted separately in
analytic_hunt.py; this FD solve is the trustworthy referee, per the overnight-run
instruction allowing qualitative/numerical verification over a shaky series.)

Nondimensional problem (duct |x|<=ar, |y|<=1; uniform field B in y; axial flow u(x,y)
driven by unit pressure gradient):

    lap(u) + Ha * d b/dy = -1        (momentum)
    lap(b) + Ha * d u/dy =  0        (induction for axial induced field b)

Boundary conditions:
    u = 0 on all four walls (no-slip)
    insulating side walls  x = +/-ar : b = 0            (no current into wall)
    conducting Hartmann walls y = +/-1 : d b/dy = 0     (perfect conductor)

Reference: Hunt, J.C.R. (1965). Magnetohydrodynamic flow in rectangular ducts.
J. Fluid Mech. 21, 577-590.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def solve_hunt(Ha, ar=1.0, nx=161, ny=161):
    """Solve the conducting-Hartmann-wall / insulating-side-wall duct.
    Returns x, y, u(ny,nx), b(ny,nx)."""
    x = np.linspace(-ar, ar, nx)
    y = np.linspace(-1.0, 1.0, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    N = nx * ny

    def idx(i, j):  # i over x (col), j over y (row)
        return j * nx + i

    # unknown vector z = [u (N), b (N)]
    A = sp.lil_matrix((2 * N, 2 * N))
    rhs = np.zeros(2 * N)

    def U(i, j): return idx(i, j)
    def B(i, j): return N + idx(i, j)

    for j in range(ny):
        for i in range(nx):
            uId, bId = U(i, j), B(i, j)
            on_side = (i == 0 or i == nx - 1)      # x = +/- ar
            on_hart = (j == 0 or j == ny - 1)      # y = +/- 1

            # ---- u equation ----
            if on_side or on_hart:
                A[uId, uId] = 1.0
                rhs[uId] = 0.0                      # no-slip u=0 on all walls
            else:
                A[uId, U(i - 1, j)] += 1.0 / dx**2
                A[uId, U(i + 1, j)] += 1.0 / dx**2
                A[uId, U(i, j - 1)] += 1.0 / dy**2
                A[uId, U(i, j + 1)] += 1.0 / dy**2
                A[uId, uId] += -2.0 / dx**2 - 2.0 / dy**2
                A[uId, B(i, j + 1)] += Ha / (2 * dy)
                A[uId, B(i, j - 1)] += -Ha / (2 * dy)
                rhs[uId] = -1.0

            # ---- b equation ----
            if on_side:
                A[bId, bId] = 1.0
                rhs[bId] = 0.0                      # insulating side walls: b=0
            elif on_hart:
                # conducting Hartmann wall: db/dy = 0 (one-sided)
                A[bId, bId] = 1.0
                if j == 0:
                    A[bId, B(i, 1)] = -1.0
                else:
                    A[bId, B(i, ny - 2)] = -1.0
                rhs[bId] = 0.0
            else:
                A[bId, B(i - 1, j)] += 1.0 / dx**2
                A[bId, B(i + 1, j)] += 1.0 / dx**2
                A[bId, B(i, j - 1)] += 1.0 / dy**2
                A[bId, B(i, j + 1)] += 1.0 / dy**2
                A[bId, bId] += -2.0 / dx**2 - 2.0 / dy**2
                A[bId, U(i, j + 1)] += Ha / (2 * dy)
                A[bId, U(i, j - 1)] += -Ha / (2 * dy)
                rhs[bId] = 0.0

    z = spla.spsolve(A.tocsr(), rhs)
    u = z[:N].reshape(ny, nx)
    b = z[N:].reshape(ny, nx)
    return x, y, u, b


def side_jet_metric(x, y, u):
    """Ratio of the near-side-wall peak velocity to the centerline velocity along
    the mid-plane y=0. >1 indicates Hunt side jets."""
    j0 = np.argmin(np.abs(y))
    line = u[j0, :]
    uc = line[np.argmin(np.abs(x))]
    return line.max() / uc, x[np.argmax(line)]


def main():
    import sys, os
    HERE = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        mpl = True
    except Exception:
        mpl = False

    print("Hunt conducting-wall duct (FD reference) — side-jet check")
    print("-" * 60)
    summary = []
    for Ha in (5, 10, 20, 50):
        x, y, u, b = solve_hunt(Ha, ar=1.0, nx=141, ny=141)
        ratio, xpk = side_jet_metric(x, y, u)
        jets = "SIDE JETS" if ratio > 1.02 else "no jets"
        print(f"  Ha={Ha:3d}: peak/centerline (y=0) = {ratio:.3f} at x={xpk:+.3f}  [{jets}]")
        summary.append((Ha, ratio, xpk))
        if mpl:
            j0 = np.argmin(np.abs(y))
            plt.figure(figsize=(6, 4))
            plt.plot(x, u[j0, :] / u[j0, np.argmin(np.abs(x))], label=f"Ha={Ha}")
            plt.xlabel("x / a (toward insulating side walls)")
            plt.ylabel("u / u_centerline at y=0")
            plt.title(f"Hunt duct mid-plane profile, Ha={Ha}")
            plt.grid(alpha=0.3); plt.legend(); plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"hunt_midplane_Ha{Ha}.png"), dpi=130)
            plt.close()
            # contour
            plt.figure(figsize=(5.2, 4.6))
            plt.contourf(x, y, u, 40)
            plt.colorbar(label="u")
            plt.xlabel("x / a"); plt.ylabel("y / a (B direction)")
            plt.title(f"Hunt duct axial velocity, Ha={Ha}")
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"hunt_contour_Ha{Ha}.png"), dpi=130)
            plt.close()

    with open(os.path.join(outdir, "hunt_sidejet_table.md"), "w") as f:
        f.write("# Hunt conducting-wall duct — side-jet metric (FD reference)\n\n")
        f.write("Peak/centerline axial velocity on the mid-plane y=0. >1 => side jets "
                "(Hunt's signature; absent in the insulating Hartmann/Shercliff case).\n\n")
        f.write("| Ha | peak/centerline | x of peak | side jets? |\n|----|----|----|----|\n")
        for Ha, r, xp in summary:
            f.write(f"| {Ha} | {r:.3f} | {xp:+.3f} | {'yes' if r>1.02 else 'no'} |\n")
    print("-" * 60)
    print("wrote", outdir)


if __name__ == "__main__":
    main()
