"""
Steady Hunt conjugate-duct x=0 centerline cross-check: MUG vs COMSOL.

Used to anchor the V2 endpoint at (Ha=100, c=1) and fill the ladder gap between
the validated Ha=20 and Ha=1323 c=1 points. Compares the Hartmann centerline cut.

MUG profile (hunt_mug.profile): header '# ...' lines then rows "x z vely by".
  core velocity nondim: u_nondim = vely * nu/(fy*a^2) = 100*vely (nu=1,fy=1e-2)
  -- BUT for the Rm-matched run fy=1e4 => u_nondim = vely*1/(1e4*1) = 1e-4*vely.
  Pass --fy to set it; u_nondim = vely*nu/(fy*a^2).
  induced field compared as by/B0 (pass --b0; MUG B0=0.1121126 at Ha=100).

COMSOL centerline CSV: '%'-comment header then rows "x,y,u,b" (x=0 cut). u is the
already-nondim axial velocity, b is by/B0 in COMSOL's B0=1 nondim.

Usage:
  python3 compare_hunt_centerline.py hunt_mug.profile \
      tier3_hunt_c1_Ha100_centerline.csv --fy 1e4 --b0 0.1121126 --nu 1
"""
import argparse

import numpy as np


def load_mug(path, fy, nu, a, b0):
    rows = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            rows.append([float(v) for v in ln.split()])
    p = np.array(rows)
    x, z, vely, by = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
    cut = np.abs(x) < 1e-6
    y = z[cut]
    u_nd = vely[cut] * nu / (fy * a ** 2)      # velocity nondim (matches COMSOL u)
    b_nd = by[cut] / b0                          # induced field as by/B0
    o = np.argsort(y)
    return y[o], u_nd[o], b_nd[o]


def load_comsol(path):
    rows = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("%"):
                continue
            parts = s.split(",")
            try:
                rows.append([float(v) for v in parts])
            except ValueError:
                continue
    p = np.array(rows)
    # columns: x, y, u, b  (x=0)
    y, u, b = p[:, 1], p[:, 2], p[:, 3]
    # COMSOL exports NaN at the no-slip walls y=+/-1 (u=0 physically) and lists
    # duplicate shared-element nodes; drop NaNs and collapse duplicate y.
    good = np.isfinite(u) & np.isfinite(b)
    y, u, b = y[good], u[good], b[good]
    yu, idx = np.unique(y, return_index=True)
    return yu, u[idx], b[idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mug")
    ap.add_argument("comsol")
    ap.add_argument("--fy", type=float, default=1e4)
    ap.add_argument("--nu", type=float, default=1.0)
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--b0", type=float, default=0.1121126)
    args = ap.parse_args()

    ym, um, bm = load_mug(args.mug, args.fy, args.nu, args.a, args.b0)
    yc, uc, bc = load_comsol(args.comsol)

    # interpolate MUG onto COMSOL's y-grid for pointwise relL2
    umi = np.interp(yc, ym, um)
    bmi = np.interp(yc, ym, bm)

    def relL2(a, b):
        return np.linalg.norm(a - b) / np.linalg.norm(b)

    # core = nearest to y=0
    j0m, j0c = np.argmin(np.abs(ym)), np.argmin(np.abs(yc))
    core_m, core_c = um[j0m], uc[j0c]
    pk_m, pk_c = np.max(np.abs(bm)), np.max(np.abs(bc))

    print("Steady Hunt c=1 centerline cross-check  (MUG vs COMSOL)")
    print(f"  core v(0,0):  MUG={core_m:.5e}  COMSOL={core_c:.5e}  "
          f"rel.err={abs(core_m-core_c)/abs(core_c):.3%}")
    print(f"  velocity relL2 v(y):  {relL2(umi, uc):.3e}")
    print(f"  peak|by/B0|:  MUG={pk_m:.4e}  COMSOL={pk_c:.4e}  "
          f"ratio={pk_m/pk_c:.4f}")
    print(f"  induced relL2 by(y):  {relL2(bmi, bc):.3e}  "
          f"(raw; if ratio!=1 the by nondim differs by that constant)")
    # also the shape-only by relL2 (normalized), Rm-independent
    bmn = bmi / (pk_m or 1); bcn = bc / (pk_c or 1)
    print(f"  induced relL2 by SHAPE (peak-normalized): {relL2(bmn, bcn):.3e}")
    verdict = ("MATCH — endpoint anchored" if abs(core_m-core_c)/abs(core_c) < 0.02
               and relL2(umi, uc) < 0.02 else "CHECK — >2%, reconcile")
    print(f"  VERDICT: {verdict}")


if __name__ == "__main__":
    main()
