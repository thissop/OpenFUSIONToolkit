"""
Gap-4 validation: LM_MHD reduced-model formulas vs. a resolved MUG duct solve.

The `OpenFUSIONToolkit.LM_MHD` module is formula-level (Hartmann drag closure,
Ohm's law, Lorentz/Joule densities, regime metrics) -- explicitly "not a solver."
This harness checks that those reduced-model formulas agree with what the resolved
MUG solver actually produces, using a resolved duct profile (x z vely by).

Two checks:
  (1) HARTMANN DRAG CLOSURE. `hartmann_channel_effective_drag` predicts the
      *insulating Hartmann channel* mean-flow response
          u_mean = f a^2/nu * [1 - tanh(Ha)/Ha] / Ha^2 .
      The resolved case is a square *duct* (side walls add drag), so the closure
      is a reduced model, not exact. We report resolved duct mean and core vs. the
      channel closure, quantifying the reduced-model error -- the honest "how good
      is the 0-D drag closure for a real duct" number. In the Ha->inf core the
      duct core approaches the channel value, so the CORE ratio is the tightest
      check.
  (2) ENERGY BALANCE. For a steady duct driven by body force f, input power
      f*<vely> must equal viscous + Ohmic dissipation. We compute all three from
      the resolved (vely, by) and check closure -- validating both the resolved
      solve's energetics and the LM_MHD power formulas.

Usage:
  python3 validate_lm_mhd_reduced.py <hunt_mug.profile> --ha 20 --fy 1e-2 --nu 1 --eta 1 --a 1
"""
import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import griddata

# import the formula-level closure directly from the module file
_cl = Path(__file__).resolve().parents[2] / "src/python/OpenFUSIONToolkit/LM_MHD/closures.py"
_spec = importlib.util.spec_from_file_location("lm_mhd_closures", _cl)
closures = importlib.util.module_from_spec(_spec)
sys.modules["lm_mhd_closures"] = closures
_spec.loader.exec_module(closures)


def load_profile(path):
    rows = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip():
                continue
            rows.append([float(v) for v in ln.split()])
    p = np.array(rows)
    return p[:, 0], p[:, 1], p[:, 2], p[:, 3]   # x, z(=y), vely, by


def to_grid(x, z, f, n=201):
    gx = np.linspace(x.min(), x.max(), n)
    gz = np.linspace(z.min(), z.max(), n)
    GX, GZ = np.meshgrid(gx, gz)
    F = griddata((x, z), f, (GX, GZ), method="linear")
    return gx, gz, F


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile")
    ap.add_argument("--ha", type=float, required=True)
    ap.add_argument("--fy", type=float, default=1e-2)
    ap.add_argument("--nu", type=float, default=1.0)
    ap.add_argument("--eta", type=float, default=1.0)
    ap.add_argument("--a", type=float, default=1.0)
    args = ap.parse_args()

    x, z, vely, by = load_profile(args.profile)
    gx, gz, V = to_grid(x, z, vely)
    _, _, B = to_grid(x, z, by)

    # --- resolved duct mean and core velocity ---
    u_mean_duct = np.nanmean(V)
    ic = np.argmin(np.abs(gx)); jc = np.argmin(np.abs(gz))
    u_core_duct = V[jc, ic]

    # --- (1) Hartmann channel drag closure ---
    alpha = closures.hartmann_channel_effective_drag(args.a, args.ha, args.nu)
    # reduced 0-D balance f = alpha * u_mean  ->  u_mean_channel = f/alpha
    u_mean_channel = args.fy / alpha
    print("=== (1) Hartmann drag closure vs resolved duct ===")
    print(f"  Ha={args.ha:g}  alpha_eff={alpha:.4e}  (nu/a^2 * Ha^2/[1-tanh Ha/Ha])")
    print(f"  channel closure   u_mean = f/alpha        = {u_mean_channel:.5e}")
    print(f"  resolved DUCT     u_mean (with side walls) = {u_mean_duct:.5e}"
          f"   ratio duct/channel = {u_mean_duct/u_mean_channel:.3f}")
    print(f"  resolved DUCT     u_core                   = {u_core_duct:.5e}"
          f"   ratio core/channel = {u_core_duct/u_mean_channel:.3f}")
    # FINDING (2026-07-14): resolved insulating DUCT core scales as 1/Ha, but the
    # channel closure scales as 1/Ha^2 -> ratio core/channel == Ha. This is real
    # physics: side walls (parallel to B) carry the return current through the side
    # layers, reducing core braking so a duct core flows ~Ha faster than a Hartmann
    # channel. The `hartmann_channel_effective_drag` closure is therefore a correct
    # CHANNEL formula but NOT a valid DUCT drag model -- it over-brakes by ~Ha.
    scaling = u_core_duct / u_mean_channel   # ~= Ha if duct core ~1/Ha
    print(f"  --> core/channel ratio = {scaling:.1f} ~= Ha={args.ha:g}: duct core "
          f"~1/Ha (Shercliff), channel ~1/Ha^2. The channel closure is NOT a duct "
          f"drag model (off by ~Ha; misses side-layer current return).")

    # --- (2) energy balance: input = viscous + Ohmic ---
    dx = gx[1]-gx[0]; dz = gz[1]-gz[0]
    dVdx = np.gradient(V, dx, axis=1); dVdz = np.gradient(V, dz, axis=0)
    dBdx = np.gradient(B, dx, axis=1); dBdz = np.gradient(B, dz, axis=0)
    m = np.isfinite(V) & np.isfinite(B)
    P_in    = args.fy * np.nanmean(V)                     # body-force input
    P_visc  = args.nu * np.nanmean(dVdx[m]**2 + dVdz[m]**2)
    # j = curl(by yhat) => |j|^2 = (dby/dx)^2+(dby/dz)^2 ; Ohmic = eta*|j|^2 (mu0=1 units)
    P_ohm   = args.eta * np.nanmean(dBdx[m]**2 + dBdz[m]**2)
    print("=== (2) steady energy balance (per unit volume, mean) ===")
    print(f"  input  f*<v>            = {P_in:.4e}")
    print(f"  viscous nu*<|grad v|^2> = {P_visc:.4e}")
    print(f"  Ohmic  eta*<|grad b|^2> = {P_ohm:.4e}")
    bal = (P_visc + P_ohm) / P_in if P_in else float('nan')
    print(f"  (viscous+Ohmic)/input   = {bal:.3f}   (1.0 => balanced)")
    print(f"  Ohmic fraction of dissipation = {P_ohm/(P_visc+P_ohm):.3f}")


if __name__ == "__main__":
    main()
