"""
Score a MUG lid-driven cavity run — the V4 feasibility gate — CORRECTLY.

History: the first version located the "vortex center" as the global interior
|v|-minimum. That is WRONG: in an under-developed cavity the global |v|-min
falls into the dead, near-stagnant LOWER cavity (a near-wall null carrying ~0
circulation), not the primary vortex. Adversarial review (2026-07-13) showed the
reported (0.49,0.12) was such an artifact while the true primary vortex sat near
Ghia at ~(0.60,0.87). This version fixes the diagnostic and adds the
incompressibility health checks the true_pressure low-Mach claim actually rests
on.

Vortex center is found as the interior extremum of the STREAMFUNCTION psi
(psi_x = -vz, psi_z = vx), cross-checked by the circulation Gamma = ∮ v·dl around
the candidate (a real vortex core has |Gamma| >> 0 of one sign). The gate also
reports, because "incompressible low-Mach" is the load-bearing claim:
  - RMS|div v| and the |div v|/|grad v| distribution (should be ~1e-3, not ~10%)
  - density deviation |n-1| if a checkpoint is supplied (low-Mach => small)
  - unsteadiness ||v(t2)-v(t1)|| if two checkpoints are supplied

Usage:
  python3 score_cavity.py cavity_mug.profile [--re 100]
          [--rst xmhd2d_04000.rst] [--rst-prev xmhd2d_03000.rst]
Profile columns: "x z velx velz". Checkpoints (HDF5) supply U_n / U_velx / U_velz
in the SAME node order as the profile rows.
"""
import argparse

import numpy as np
from scipy.interpolate import griddata

# Ghia et al. (1982) primary-vortex center, unit cavity, lid on top
GHIA = {100: (0.6172, 0.7344), 400: (0.5547, 0.6055), 1000: (0.5313, 0.5625)}


def grid_fields(x, z, vx, vz, n=160, lo=0.05, hi=0.95):
    gx = np.linspace(lo, hi, n)
    gz = np.linspace(lo, hi, n)
    GX, GZ = np.meshgrid(gx, gz)
    VX = griddata((x, z), vx, (GX, GZ), method="linear")
    VZ = griddata((x, z), vz, (GX, GZ), method="linear")
    return gx, gz, GX, GZ, VX, VZ


def streamfunction(gx, gz, VX, VZ):
    """psi with d(psi)/dz = vx, d(psi)/dx = -vz. Integrate vx along z (per column)
    from the bottom; robust enough to locate the interior extremum."""
    dz = gz[1] - gz[0]
    psi = np.zeros_like(VX)
    # cumulative integral of vx up each column (axis 0 is z)
    psi = np.cumsum(np.nan_to_num(VX), axis=0) * dz
    psi[~np.isfinite(VX)] = np.nan
    return psi


def circulation(gx, gz, GX, GZ, VX, VZ, cx, cz, r=0.15):
    """Gamma = ∮ v·dl on a circle radius r about (cx,cz), sampled from the grid."""
    th = np.linspace(0, 2 * np.pi, 128, endpoint=False)
    px, pz = cx + r * np.cos(th), cz + r * np.sin(th)
    vx = griddata((GX.ravel(), GZ.ravel()), VX.ravel(), (px, pz), method="linear")
    vz = griddata((GX.ravel(), GZ.ravel()), VZ.ravel(), (px, pz), method="linear")
    # tangential component * dl ; dl = r dth, tangent = (-sin, cos)
    tang = -vx * np.sin(th) + vz * np.cos(th)
    return np.nansum(tang) * r * (th[1] - th[0])


def divergence_health(gx, gz, VX, VZ):
    dx, dz = gx[1] - gx[0], gz[1] - gz[0]
    dvx_dx = np.gradient(VX, dx, axis=1)
    dvz_dz = np.gradient(VZ, dz, axis=0)
    div = dvx_dx + dvz_dz
    dvx_dz = np.gradient(VX, dz, axis=0)
    dvz_dx = np.gradient(VZ, dx, axis=1)
    grad = np.sqrt(dvx_dx ** 2 + dvz_dz ** 2 + dvx_dz ** 2 + dvz_dx ** 2)
    m = np.isfinite(div) & np.isfinite(grad) & (grad > 0)
    ratio = np.abs(div[m]) / grad[m]
    return np.sqrt(np.nanmean(div[m] ** 2)), np.percentile(ratio, [50, 95])


def load_rst(path, key):
    import h5py
    with h5py.File(path, "r") as f:
        return np.asarray(f[key]).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--re", type=int, default=100)
    ap.add_argument("--rst", help="checkpoint for density / incompressibility")
    ap.add_argument("--rst-prev", help="earlier checkpoint for unsteadiness")
    args = ap.parse_args()

    pts = np.loadtxt(args.profile)
    x, z, vx, vz = pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3]
    ulid = np.nanmax(np.abs(vx))
    gx, gz, GX, GZ, VX, VZ = grid_fields(x, z, vx, vz)

    # --- vortex center via streamfunction extremum (interior only) ---
    psi = streamfunction(gx, gz, VX, VZ)
    interior = (GX > 0.1) & (GX < 0.9) & (GZ > 0.1) & (GZ < 0.9)
    pmask = np.where(interior & np.isfinite(psi), np.abs(psi), -np.inf)
    j = np.unravel_index(np.argmax(pmask), psi.shape)
    cx, cz = GX[j], GZ[j]
    Gam = circulation(gx, gz, GX, GZ, VX, VZ, cx, cz)

    # recirculation: is there return flow (vx<0) in the interior?
    has_return = np.nanmin(VX[interior]) < -0.02 * ulid
    rms_div, div_ratio = divergence_health(gx, gz, VX, VZ)

    print(f"Lid-driven cavity, Re={args.re}  (u_lid={ulid:.3e}, {pts.shape[0]} pts)")
    print(f"  return flow present (vx<0 interior): {has_return}")
    print(f"  vortex center (streamfunction extremum): x={cx:.3f}, z={cz:.3f}")
    print(f"  circulation about center (r=0.15): {Gam:.3e}  "
          f"(|Gamma|>>0 one-sign => real vortex core)")
    if args.re in GHIA:
        gxc, gzc = GHIA[args.re]
        print(f"  Ghia Re={args.re} center: x={gxc:.3f}, y={gzc:.3f}  "
              f"-> distance {np.hypot(cx - gxc, cz - gzc):.3f}")

    print(f"  INCOMPRESSIBILITY: RMS|div v| = {rms_div:.2e}  "
          f"({rms_div/ (ulid):.2e} in u_lid/L)")
    print(f"    |div v|/|grad v| median={div_ratio[0]:.3f}, 95pct={div_ratio[1]:.3f}"
          f"   (clean incompressible ~<0.01)")

    n_ok = True
    if args.rst:
        n = load_rst(args.rst, "U_n")
        ndev = np.max(np.abs(n - 1.0))
        print(f"  density deviation max|n-1| = {ndev:.3f}  "
              f"(low-Mach M~0.2 expects ~0.03)")
        n_ok = ndev < 0.05
        if args.rst_prev:
            npv = load_rst(args.rst_prev, "U_n")
            print(f"    n-std growth: {np.std(npv-1):.2e} -> {np.std(n-1):.2e}")
    if args.rst and args.rst_prev:
        vxp = load_rst(args.rst_prev, "U_velx"); vzp = load_rst(args.rst_prev, "U_velz")
        vxc = load_rst(args.rst, "U_velx"); vzc = load_rst(args.rst, "U_velz")
        sp = np.hypot(vxp, vzp); sc = np.hypot(vxc, vzc)
        unsteady = np.linalg.norm(sc - sp) / np.linalg.norm(sc)
        print(f"  UNSTEADINESS (speed field, consecutive rst): {unsteady:.3f}  "
              f"(steady => ~0)")

    # honest verdict
    recirc = has_return and abs(Gam) > 0.02 * ulid
    incompressible = (div_ratio[0] < 0.03) and n_ok
    if recirc and incompressible:
        v = "V4 FEASIBLE — clean incompressible recirculation"
    elif recirc:
        v = ("V4 QUALIFIED — a recirculating cell forms, but low-Mach "
             "incompressibility NOT clean (see div / n) => held")
    else:
        v = "V4 NO-GO — no coherent recirculation"
    print(f"  VERDICT: {v}")


if __name__ == "__main__":
    main()
