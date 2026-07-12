"""
Score the MUG V3 entrance-flow run against the COMSOL reference fields
(lm-mhd-comsol data/results/v3_entrance_Ha100.csv + v3_centerline.csv).

Mapping: MUG in-plane channel (x streamwise, z wall-normal, flow velx, flux
function psi) -> COMSOL (X streamwise, Y wall-normal, u, b). The streamwise
induced field is b = +/- d(psi)/dz; the sign is fixed empirically by matching
COMSOL's b (both codes' b enter quadratically in the physics, so the relative
overlay is what matters).

Nondim (f = +1 units, fx = 1e-2, nu = eta = 1, a = 1):
  u* = velx / fx
  b* = (d psi / dz) / S,  S = fx * a^2 * sqrt(mu0*rho) / sqrt(eta*nu)

Diagnostics: field relL2(u), relL2(b), centerline overlay, L_e(99%) both codes.

Usage: python3 compare_v3.py v3_entrance_mug.profile v3_entrance_Ha100.csv v3_centerline.csv
"""
import sys

import numpy as np
from scipy.interpolate import griddata

MU0 = 4e-7 * np.pi
M_P = 1.67262192369e-27


def main():
    mug_path, field_csv, cl_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    fx = 1e-2
    rho = M_P * 5.98e26
    S = fx * np.sqrt(MU0 * rho)

    pts = np.loadtxt(mug_path)
    x, z, u, psi = pts[:, 0], pts[:, 1], pts[:, 2] / fx, pts[:, 3]

    # psi -> b on a regular grid (interior; avoid wall extrapolation)
    xg = np.linspace(0.01, 3.99, 400)
    zg = np.linspace(-0.995, 0.995, 240)
    XX, ZZ = np.meshgrid(xg, zg)
    PSI = griddata((x, z), psi, (XX, ZZ), method="linear")
    UU = griddata((x, z), u, (XX, ZZ), method="linear")
    Bg = np.gradient(PSI, zg, axis=0) / S

    ref = np.genfromtxt(field_csv, delimiter=",", comments="%")
    good = np.isfinite(ref[:, 2])
    Xr, Yr, Ur, Br = ref[good, 0], ref[good, 1], ref[good, 2], ref[good, 3]
    ui = griddata((XX.ravel(), ZZ.ravel()), UU.ravel(),
                  np.column_stack([Xr, Yr]), method="linear")
    bi = griddata((XX.ravel(), ZZ.ravel()), Bg.ravel(),
                  np.column_stack([Xr, Yr]), method="linear")
    ok = np.isfinite(ui) & np.isfinite(bi)
    relu = np.sqrt(np.nansum((ui[ok] - Ur[ok]) ** 2) / np.nansum(Ur[ok] ** 2))
    # b sign: pick the sign minimizing the residual (flux convention freedom)
    r_plus = np.nansum((bi[ok] - Br[ok]) ** 2)
    r_minus = np.nansum((-bi[ok] - Br[ok]) ** 2)
    sgn = 1.0 if r_plus <= r_minus else -1.0
    relb = np.sqrt(min(r_plus, r_minus) / np.nansum(Br[ok] ** 2))

    cl = np.genfromtxt(cl_csv, delimiter=",", comments="%")
    xc, uc_ref = cl[:, 0], cl[:, 2]
    uc_mug = griddata((x, z), u, np.column_stack([xc, np.zeros_like(xc)]),
                      method="linear")

    def le99(xv, uv):
        u_inf = uv[-1]
        above = np.where(uv >= 0.99 * u_inf)[0]
        return xv[above[0]] if above.size else np.nan

    print(f"V3 entrance, Ha=100 ({pts.shape[0]} MUG pts, {ok.sum()} ref pts):")
    print(f"  field relL2(u) = {relu:.3e}")
    print(f"  field relL2(b) = {relb:.3e}   (b sign = {sgn:+.0f} d(psi)/dz)")
    print(f"  u_c(L): MUG = {uc_mug[-1]:.6e}, COMSOL = {uc_ref[-1]:.6e}")
    print(f"  L_e(99%): MUG = {le99(xc, uc_mug):.3f}, COMSOL = {le99(xc, uc_ref):.3f}")
    print(f"  centerline max rel dev = "
          f"{np.nanmax(np.abs(uc_mug - uc_ref) / np.nanmax(uc_ref)):.3e}")


if __name__ == "__main__":
    main()
