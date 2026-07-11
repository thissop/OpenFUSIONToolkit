"""
Score the MUG V4a trapezoid-duct run against the COMSOL reference
(lm-mhd-comsol data/results/v4a_trapduct_Ha100.csv; targets Q = 2.76810e-2,
max|b| = 9.00e-3 in forcing=-1 nondim units).

Nondimensionalization (matches SYNC contract):
  u* = vely * nu / (fy * a^2)
  b* = by / S,  S = fy * a^2 * sqrt(mu0*rho) / sqrt(eta*nu)
(derived by matching MUG's dimensional induction balance to COMSOL's
 lap(b*) + Ha * du*/dy = 0 system; rho = m_p * n0 * den_scale.)

Flow rate: Q = int u* dA over the trapezoid via Delaunay triangulation of the
mesh points (boundary u=0 included, so the integral is well-posed).

Usage: python3 compare_trapduct.py trapduct_mug.profile v4a_trapduct_Ha100.csv
"""
import re
import sys

import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import Delaunay

MU0 = 4e-7 * np.pi
M_P = 1.67262192369e-27


def load_mug(path, nu=1.0, eta=1.0, a=1.0, n0=1.0, den_scale=5.98e26):
    Ha = fy = None
    rows = []
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
            c = line.split()
            if len(c) == 4:
                rows.append([float(v) for v in c])
    pts = np.array(rows)
    rho = M_P * n0 * den_scale
    S = fy * a**2 * np.sqrt(MU0 * rho) / np.sqrt(eta * nu)
    u = pts[:, 2] * nu / (fy * a**2)
    b = pts[:, 3] / S
    return Ha, pts[:, 0], pts[:, 1], u, b


def flow_rate(x, z, u):
    tri = Delaunay(np.column_stack([x, z]))
    s = tri.simplices
    p = tri.points
    v1 = p[s[:, 1]] - p[s[:, 0]]
    v2 = p[s[:, 2]] - p[s[:, 0]]
    area = 0.5 * np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0])
    return float(np.sum(area * u[s].mean(axis=1)))


def main():
    mug_path, csv_path = sys.argv[1], sys.argv[2]
    Ha, x, z, u, b = load_mug(mug_path)
    Q = flow_rate(x, z, u)
    u0 = float(griddata((x, z), u, np.array([[0.0, 0.0]]), method="linear")[0])
    print(f"MUG V4a (Ha={Ha:.4g}, {x.size} points):")
    print(f"  Q  = {Q:.6e}   (COMSOL target 2.76810e-2, err "
          f"{abs(Q - 2.76810e-2)/2.76810e-2*100:.3f}%)")
    print(f"  u(0,0) = {u0:.6e},  Ha*u0 = {Ha*u0:.5f}")
    print(f"  max|b| = {np.max(np.abs(b)):.4e}   (COMSOL target 9.00e-3)")

    # field-level comparison on the COMSOL grid (NaN marks outside-trapezoid)
    ref = np.genfromtxt(csv_path, delimiter=",", comments="%")
    good = np.isfinite(ref[:, 2])
    Xr, Yr, Ur, Br = ref[good, 0], ref[good, 1], ref[good, 2], ref[good, 3]
    ui = griddata((x, z), u, np.column_stack([Xr, Yr]), method="linear")
    bi = griddata((x, z), b, np.column_stack([Xr, Yr]), method="linear")
    ok = np.isfinite(ui)
    relu = np.sqrt(np.nansum((ui[ok] - Ur[ok]) ** 2) / np.nansum(Ur[ok] ** 2))
    relb = np.sqrt(np.nansum((bi[ok] - Br[ok]) ** 2) / np.nansum(Br[ok] ** 2))
    print(f"  field relL2(u) vs COMSOL = {relu:.3e}  ({ok.sum()} ref points)")
    print(f"  field relL2(b) vs COMSOL = {relb:.3e}")


if __name__ == "__main__":
    main()
