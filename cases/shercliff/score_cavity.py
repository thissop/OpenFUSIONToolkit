"""
Score a MUG lid-driven cavity run — the V4 feasibility gate.

Question: does MUG's true-pressure low-Mach mode produce a recirculating primary
vortex (the capability the V4 backward-facing step needs)? If yes, and the vortex
center is near the Ghia et al. (1982) Re=100 benchmark (~x=0.617, y=0.738 in a
unit cavity), V4 is feasible. If the interior flow doesn't recirculate (no
stagnation point away from the lid), V4 is off on MUG.

Usage: python3 score_cavity.py cavity_mug.profile [--re 100]
"""
import argparse

import numpy as np
from scipy.interpolate import griddata

# Ghia et al. (1982) primary-vortex center, unit cavity, lid on top
GHIA = {100: (0.6172, 0.7344), 400: (0.5547, 0.6055), 1000: (0.5313, 0.5625)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--re", type=int, default=100)
    args = ap.parse_args()

    pts = np.loadtxt(args.profile)
    x, z, vx, vz = pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3]
    speed = np.hypot(vx, vz)
    ulid = np.nanmax(np.abs(vx))

    # interior grid (exclude walls and the lid boundary layer)
    gx = np.linspace(0.15, 0.9, 160)
    gz = np.linspace(0.1, 0.85, 160)
    GX, GZ = np.meshgrid(gx, gz)
    S = griddata((x, z), speed, (GX, GZ), method="linear")
    VX = griddata((x, z), vx, (GX, GZ), method="linear")

    # recirculation signature: is there return flow (vx<0) in the lower interior?
    lower = GZ < 0.5
    has_return = np.nanmin(VX[lower]) < -0.02 * ulid
    # vortex center = interior speed minimum (stagnation of the recirculation)
    j = np.nanargmin(np.where(np.isfinite(S), S, np.inf))
    cx, cz = GX.ravel()[j], GZ.ravel()[j]
    smin = S.ravel()[j]

    print(f"Lid-driven cavity, Re={args.re}  (u_lid={ulid:.3e}, {pts.shape[0]} pts)")
    print(f"  return flow present (vx<0 in lower half): {has_return}")
    print(f"  min interior |v|/u_lid = {smin/ulid:.3e}  (stagnation => vortex core)")
    print(f"  vortex center (MUG): x={cx:.3f}, z={cz:.3f}")
    if args.re in GHIA:
        gxc, gzc = GHIA[args.re]
        d = np.hypot(cx - gxc, cz - gzc)
        print(f"  Ghia Re={args.re} center: x={gxc:.3f}, y={gzc:.3f}  -> distance {d:.3f}")
    verdict = ("V4 FEASIBLE — recirculation present" if has_return and smin < 0.1 * ulid
               else "V4 NO-GO — no clear recirculation")
    print(f"  VERDICT: {verdict}")


if __name__ == "__main__":
    main()
