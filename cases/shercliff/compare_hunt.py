"""
Score a MUG Hunt-duct run (test_hunt_mug.F90) against the FD referee
(cases/hunt/hunt_duct.py — perfect-conductor Hartmann walls, insulating sides).

Mapping: MUG mesh (x = side direction, z = Hartmann direction, B along z)
-> referee (x = side |x|<=ar insulating, y = Hartmann |y|<=1 conducting).
Nondim: u* = vely*nu/(fy*a^2) (unit-forcing convention, same as Shercliff).

Signature physics vs Shercliff: high-velocity side-wall jets — report the
jet peak ratio u_max/u(0,0) as well as core and relL2.

Usage: python3 compare_hunt.py hunt_mug.profile [--nfd 241]
"""
import argparse
import re
import sys

import numpy as np
from scipy.interpolate import griddata

sys.path.insert(0, __file__.rsplit("/", 2)[0] + "/hunt")
from hunt_duct import solve_hunt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--nfd", type=int, default=241)
    ap.add_argument("--nu", type=float, default=1.0)
    ap.add_argument("--a", type=float, default=1.0)
    args = ap.parse_args()

    Ha = fy = None
    rows = []
    with open(args.profile) as fh:
        for line in fh:
            if line.startswith("#"):
                m = re.search(r"Ha_pred\s*=\s*([0-9.eE+-]+)", line)
                if m:
                    Ha = float(m.group(1))
                m = re.search(r"fy\s*=\s*([0-9.eE+-]+)", line)
                if m:
                    fy = float(m.group(1))
                continue
            c = line.split()
            if len(c) >= 3:
                rows.append([float(v) for v in c[:3]])
    pts = np.array(rows)
    u_mug = pts[:, 2] * args.nu / (fy * args.a**2)

    xf, yf, u_fd, _ = solve_hunt(Ha, ar=1.0, nx=args.nfd, ny=args.nfd)
    XX, YY = np.meshgrid(xf, yf)
    keep = (np.abs(XX) <= 0.98) & (np.abs(YY) <= 0.98)
    ui = griddata((pts[:, 0], pts[:, 1]), u_mug,
                  np.column_stack([XX[keep], YY[keep]]), method="linear")
    ok = np.isfinite(ui)
    ufk = u_fd[keep][ok]
    relL2 = np.sqrt(np.sum((ui[ok] - ufk) ** 2) / np.sum(ufk ** 2))

    u0_mug = float(griddata((pts[:, 0], pts[:, 1]), u_mug,
                            np.array([[0.0, 0.0]]), method="linear")[0])
    u0_fd = float(u_fd[args.nfd // 2, args.nfd // 2])
    jet_mug = float(np.max(ui[ok])) / u0_mug
    jet_fd = float(np.max(ufk)) / u0_fd

    print(f"Hunt duct, Ha = {Ha:g} ({pts.shape[0]} MUG points, "
          f"FD referee {args.nfd}x{args.nfd}):")
    print(f"  u(0,0): MUG = {u0_mug:.6e}, FD = {u0_fd:.6e}, "
          f"err = {abs(u0_mug-u0_fd)/u0_fd*100:.3f}%")
    print(f"  relL2(u) interior = {relL2:.3e}")
    print(f"  side-jet peak / core: MUG = {jet_mug:.4f}, FD = {jet_fd:.4f}")


if __name__ == "__main__":
    main()
