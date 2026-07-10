"""
Score a MUG Shercliff duct run against the exact series (VALIDATION_CASES.md S1/S1a).

Reads shercliff_mug.profile ("x z vely" over mesh vertices, headers carry Ha_pred
and fy), nondimensionalizes u = vely*nu/(fy*a^2) (unit-forcing convention agreed
in SYNC/comsol_to_oft.md Q1), and reports the contract metrics: core-velocity
error, relL2(u) over the cross-section, and Ha*U0.

Usage: python3 compare_mug.py path/to/shercliff_mug.profile [--nu 1.0] [--a 1.0]
"""
import argparse
import re

import numpy as np
from scipy.interpolate import griddata

from shercliff_series import solve_series, core_velocity


def load_profile(path):
    Ha = fy = None
    pts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                m = re.search(r"Ha_pred\s*=\s*([0-9.eE+-]+)", line)
                if m:
                    Ha = float(m.group(1))
                m = re.search(r"fy\s*=\s*([0-9.eE+-]+)", line)
                if m:
                    fy = float(m.group(1))
                continue
            cols = line.split()
            if len(cols) == 3:
                pts.append([float(c) for c in cols])
    if Ha is None or fy is None:
        raise ValueError(f"{path}: missing Ha_pred/fy header")
    return Ha, fy, np.array(pts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--nu", type=float, default=1.0)
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--nmodes", type=int, default=4000)
    args = ap.parse_args()

    Ha, fy, pts = load_profile(args.profile)
    x_m, z_m, vely = pts[:, 0], pts[:, 1], pts[:, 2]
    u_mug = vely * args.nu / (fy * args.a**2)

    # analytic on a grid clear of the walls; series coords: x = side dir, y = Hartmann dir.
    # MUG mapping: mesh x = side direction, mesh z = Hartmann direction (B_0 along z).
    xg = np.linspace(-0.98, 0.98, 197)
    yg = np.linspace(-0.98, 0.98, 197)
    _, _, u_an, _ = solve_series(Ha, 1.0, args.nmodes, x=xg, y=yg)
    XX, YY = np.meshgrid(xg, yg)
    u_i = griddata((x_m, z_m), u_mug, (XX, YY), method="linear")
    ok = np.isfinite(u_i)

    relL2 = np.sqrt(np.sum((u_i[ok] - u_an[ok]) ** 2) / np.sum(u_an[ok] ** 2))
    U0_an = core_velocity(Ha, 1.0, args.nmodes)
    U0_mug = float(griddata((x_m, z_m), u_mug, np.array([[0.0, 0.0]]),
                            method="linear")[0])

    print(f"Ha = {Ha:.6g}   (fy = {fy:g}, {pts.shape[0]} mesh points)")
    print(f"  U0  analytic = {U0_an:.8e}")
    print(f"  U0  MUG      = {U0_mug:.8e}")
    print(f"  core-vel err = {abs(U0_mug - U0_an) / U0_an * 100:.3f}%")
    print(f"  relL2(u)     = {relL2:.3e}   (interior |x|,|z| <= 0.98)")
    print(f"  Ha * U0_mug  = {Ha * U0_mug:.6f}")


if __name__ == "__main__":
    main()
